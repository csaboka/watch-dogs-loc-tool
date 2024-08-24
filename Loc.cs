using System;
using System.Collections.Generic;
using System.IO;
using System.IO.MemoryMappedFiles;
using System.Linq;
using System.Text;
using Gibbed.IO;

namespace watch_dogs_loc
{
    public class Loc
    {
        short magic, version;
        short language;
        ushort table_length;
        uint tree_offset;
        Table[] table;

        uint tree_meta_offset;

        uint[] tree_meta;

        uint nodes_count;

        uint[] tree_entries;

        public void Read(Stream input, MemoryMappedViewAccessor accessor)
        {
            magic = input.ReadValueS16();
            version = input.ReadValueS16();
            if (magic != 0x4C53 /*LS*/ || version != 1)
            {
                Console.WriteLine("Not a valid loc file!");
                Environment.Exit(1);
            }
            language = input.ReadValueS16();
            table_length = input.ReadValueU16();
            tree_offset = input.ReadValueU32();
            nodes_count = accessor.ReadUInt32(tree_offset);
            tree_entries = new uint[nodes_count];
            accessor.ReadArray(tree_offset, tree_entries, 0, (int)nodes_count);
            List<uint> tree_meta_tmp = new List<uint>();
            uint ptr = tree_offset;
            do
            {
                tree_meta_tmp.Add(accessor.ReadUInt32(ptr - 4));
                tree_meta_tmp.Add(accessor.ReadUInt32(ptr - 8));
                ptr -= 8;
            } while (tree_meta_tmp[tree_meta_tmp.Count - 2] != 0xFFFFFFFF);
            tree_meta = tree_meta_tmp.ToArray();
            table = new Table[table_length];
            long end_of_tables = -1;
            for (ushort i = 0; i < table_length; ++i)
            {
                table[i] = new Table();
                end_of_tables = table[i].Read(input, accessor, this);
            }
            uint end_of_tables_aligned = (uint)(end_of_tables + (-end_of_tables & 3));
            if (end_of_tables_aligned != ptr)
            {
                Console.WriteLine("Unexpected data at " + end_of_tables_aligned);
                Environment.Exit(1);
            }
        }

        public void Export(String filename)
        {
            using (StreamWriter text = new StreamWriter(filename + ".txt", false, System.Text.Encoding.Unicode))
            {
                foreach (Table t in table)
                {
                    foreach (SubTableIds ids in t.sub_table_ids)
                    {
                        foreach (List<Id> ids_64 in ids.ids)
                        {
                            foreach (Id id in ids_64)
                            {
                                if (id.is_pseudo)
                                {
                                    continue;
                                }
                                string line = DecodeString(id.tree_pointers);
                                id.str = line;
                                text.WriteLine(id.id + "=" + line.Replace("\r", "[CR]").Replace("\n", "[LF]"));
                            }
                        }
                    }
                }
            }
        }

        private string DecodeString(uint[] tree_positions)
        {
            Stack<uint> stack = new Stack<uint>(tree_positions.Reverse());
            StringBuilder builder = new StringBuilder();
            while (stack.Count > 0)
            {
                uint idx = stack.Pop();
                uint tree_node = tree_entries[idx];
                if (tree_node <= 0xFFFF)
                {
                    builder.Append((char)tree_node);
                }
                else
                {
                    stack.Push(tree_node & 0xFFFF);
                    stack.Push(tree_node >> 16);
                }
            }
            return builder.ToString();
        }

        public uint FindTreePosition(BitReader reader, ref int bit_length)
        {
            uint current_uint_masked = reader.PeekBits(27);
            int bits_to_read;
            uint offset;
            int i = 0;
            while (current_uint_masked >= tree_meta[i])
            {
                i += 2;
            }
            offset = tree_meta[i + 1];
            bits_to_read = (int)(offset & 0x1F);
            bit_length -= bits_to_read;
            return (reader.ReadBits(bits_to_read) + offset) >> (32 - bits_to_read);
        }

    }

    public class Table
    {
        public uint first_id;
        public uint offset_length; // 28 bits offset + 4 bits length

        public uint offset;
        public uint length;

        public SubTableMeta[] sub_table_metas;
        public SubTableIds[] sub_table_ids;

        public long Read(Stream input, MemoryMappedViewAccessor accessor, Loc loc)
        {
            first_id = input.ReadValueU32();
            offset_length = input.ReadValueU32();

            long save = input.Position;

            offset = offset_length >> 4;
            length = offset_length & 15;

            // read meta data for each subtable
            input.Position = offset;
            sub_table_metas = new SubTableMeta[length];
            for (uint i = 0; i < length; ++i)
            {
                sub_table_metas[i] = new SubTableMeta();
                sub_table_metas[i].Read(input);
            }
            // read ids and bit offset (lo and hi) for each subtable
            sub_table_ids = new SubTableIds[length];
            uint block_first_id = first_id;
            long block_position = input.Position;
            for (uint i = 0; i < length; ++i)
            {
                input.Position = block_position;

                block_first_id += sub_table_metas[i].delta_from_prev_id;

                sub_table_ids[i] = new SubTableIds((uint)block_position);
                sub_table_ids[i].Read(ref block_first_id, sub_table_metas[i], input, accessor, loc);

                block_position += sub_table_metas[i].size;
            }

            input.Position = save;

            return block_position;
        }
    }

    public class SubTableMeta
    {
        public uint max_id, size, delta_from_prev_id;

        public void Read(Stream input)
        {
            uint first = input.ReadValueU16();
            uint second = input.ReadValueU16();

            uint whole = (first << 16) + second;
            delta_from_prev_id = second;
            if (whole >= 0x80000000)
            {
                size = input.ReadValueU16();
                if (((whole >> 30) & 1) != 0)
                {
                    uint extra = input.ReadValueU16();
                    delta_from_prev_id += (extra << 16);
                }
                max_id = first & 0x3FFF;
            }
            else
            {
                max_id = first >> 7;
                size = (whole >> 12) & 0x7FF;
                delta_from_prev_id &= 0xFFF;
            }
        }
    }

    public class SubTableIds
    {
        public List<List<Id>> ids;

        public uint start;

        public byte[] raw;

        public SubTableIds(uint start)
        {
            this.start = start;
        }

        public void Read(ref uint id_begin, SubTableMeta subTable, Stream input, MemoryMappedViewAccessor accessor, Loc loc)
        {
            ids = new List<List<Id>>();

            long subtable_ids_begin = input.Position;

            uint id_count = subTable.max_id + 1;
            uint block_64ids_count = (id_count - 1) >> 6;

            ushort[]  block_64ids_offsets = new ushort[block_64ids_count];
            for (uint i = 0; i < block_64ids_count; ++i)
            {
                block_64ids_offsets[i] = input.ReadValueU16();
            }

            for (uint j = 0; j < id_count; j += 64)
            {
                List<Id> ids_in_block = new List<Id>();
                if (j >= 64 && j % 64 == 0)
                {
                    long block_64ids_offset = block_64ids_offsets[(j >> 6) - 1];
                    if (block_64ids_offset == 0)
                    {
                        ids.Add(ids_in_block);
                        continue;
                    }
                    input.Position = subtable_ids_begin + block_64ids_offset;
                }

               
                uint current_size_in_bits = 0;
                uint k = 0;
                while (k < Math.Min(id_count - j, 64))
                {
                    Id new_id = new Id(id_begin + j + k);
                    new_id.Read(ref k, ref current_size_in_bits, input);
                    ids_in_block.Add(new_id);
                }
                BitReader bitReader = new BitReader(accessor, (int)(input.Position * 8));
                foreach (Id id in ids_in_block)
                {
                    if (id.is_pseudo)
                    {
                        continue;
                    }
                    List<uint> tree_pointers = new List<uint>();
                    int remaining_bits = (int)id.increment;
                    while (remaining_bits > 0)
                    {
                        uint entry = loc.FindTreePosition(bitReader, ref remaining_bits);
                        tree_pointers.Add(entry);
                    }
                    bitReader.Seek(remaining_bits);
                    id.tree_pointers = tree_pointers.ToArray();
                }
                ids.Add(ids_in_block);
            }

            id_begin += id_count;
        }
    }

    public class Id
    {
        // if is_pseudo then id increment
        // else size increment in bits (encoded.Length << 3)
        public uint increment;

        // calculated
        public uint id;


        public String str;
        public uint[] tree_pointers;
        public bool is_pseudo;

        public Id(uint id)
        {
            this.id = id;
            this.is_pseudo = false;
            this.str = null;
        }

        public void Read(ref uint k, ref uint current_size_in_bits, Stream input)
        {
            uint current_size = input.ReadValueU8();
            if (current_size > 0xF0)
            {
                increment = current_size - 240;
                k += (current_size - 240);
                if (k > 64)
                {
                    Console.WriteLine("Extras! last id=" + id + " skip=" + (current_size - 240) + " extra=" + (k - 64));
                    Environment.Exit(1);
                }
                is_pseudo = true;
            }
            else
            {
                if (current_size == 0xF0)
                {
                    current_size = input.ReadValueU8();
                    current_size = ((current_size << 8) + input.ReadValueU8()) + 5340;
                }
                else if (current_size >= 0xDC)
                {
                    current_size = ((current_size << 8) + input.ReadValueU8()) - 56100;
                }
                if (current_size != 0)
                {
                    current_size = 2 * current_size + 4;
                }
                increment = current_size;
                current_size_in_bits += current_size;
                ++k;
            }
        }
    }

    public class BitReader
    {
        readonly MemoryMappedViewAccessor accessor;
        int position;

        public BitReader(MemoryMappedViewAccessor accessor, int position)
        {
            this.accessor = accessor;
            this.position = position;
        }

        public uint PeekBits(int count)
        {
            if (count < 1 || count > 32)
            {
                throw new ArgumentException("Bit count out of range: " + count);
            }
            ulong answer = accessor.ReadUInt64(position >> 3).Swap();
            answer >>= (32 - (position & 7));
            int mask = int.MinValue >> (count - 1);
            return ((uint)answer) & ((uint)mask);
        }

        public uint ReadBits(int count)
        {
            uint answer = PeekBits(count);
            this.position += count;
            return answer;
        }

        public void Seek(int offset)
        {
            this.position += offset;
        }
    }
}
