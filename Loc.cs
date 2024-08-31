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
        short language;
        Table[] table;
        uint[] tree_meta;
        uint[] tree_entries;

        public void Read(Stream input, MemoryMappedViewAccessor accessor)
        {
            short magic = input.ReadValueS16();
            short version = input.ReadValueS16();
            if (magic != 0x4C53 /*LS*/ || version != 1)
            {
                Console.WriteLine("Not a valid loc file!");
                Environment.Exit(1);
            }
            language = input.ReadValueS16();
            ushort table_length = input.ReadValueU16();
            uint tree_offset = input.ReadValueU32();
            uint nodes_count = accessor.ReadUInt32(tree_offset);
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

        public void Write(Stream output)
        {
            output.WriteValueS16(0x4c53);
            output.WriteValueS16(1);
            output.WriteValueS16(language);

            List<Table> newTables = new List<Table>();
            MemoryStream tableData = new MemoryStream();
            foreach (Table oldTable in table)
            {
                List<Table> newlyAdded = oldTable.WriteData(tableData, this);
                newTables.AddRange(newlyAdded);
            }
            output.WriteValueS16((short)newTables.Count);
            long tree_offset_pos = output.Position;
            output.WriteValueU32(0);
            uint adjust = (uint)(output.Position + 8 * newTables.Count);
            foreach (Table table in newTables)
            {
                table.offset += adjust;
                table.Write(output);
            }
            output.WriteBytes(tableData.ToArray());
            while ((output.Position & 3) != 0)
            {
                output.WriteByte(0);
            }
            foreach (uint meta in tree_meta.Reverse())
            {
                output.WriteValueU32(meta);
            }
            long tree_offset = output.Position;
            foreach (uint entry in tree_entries)
            {
                output.WriteValueU32(entry);
            }
            output.Position = tree_offset_pos;
            output.WriteValueU32((uint)tree_offset);
        }

        public void Export(String filename)
        {
            using (StreamWriter text = new StreamWriter(filename + ".txt", false, System.Text.Encoding.Unicode))
            {
                foreach (Id id in AllIds())
                {
                    string line = DecodeString(id.tree_pointers);
                    id.str = line;
                    text.WriteLine(id.id + "=" + line.Replace("\r", "[CR]").Replace("\n", "[LF]"));
                }
            }
        }

        private IEnumerable<Id> AllIds()
        {

            foreach (Table t in table)
            {
                foreach (SubTableIds ids in t.sub_table_ids)
                {
                    foreach (List<Id> ids_64 in ids.ids)
                    {
                        foreach (Id id in ids_64)
                        {
                            if (!id.is_pseudo)
                            {
                                yield return id;
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

        public void WriteTreePosition(BitWriter writer, uint value)
        {
            for (int i = 0; ; i += 2)
            {
                uint offset = tree_meta[i + 1];
                int bits_to_write = (int)(offset & 0x1F);
                uint candidate = (value << (32 - bits_to_write)) - (offset & ~0x1fu);
                bool matchesEarlier = false;
                for (int j = 0; j < i; j += 2)
                {
                    matchesEarlier |= candidate < tree_meta[j];
                }
                if (!matchesEarlier && candidate < tree_meta[i] && ((candidate + offset) >> (32 - bits_to_write)) == value)
                {
                    // found the right offset and bit count
                    writer.Write(candidate >> (32 - bits_to_write), bits_to_write);
                    return;
                }
            }
        }

        public void Update(Dictionary<uint, string> newStrings)
        {

            foreach (Id id in AllIds())
            {
                string newValue;
                if (newStrings.TryGetValue(id.id, out newValue))
                {
                    id.str = newValue;
                }
            }
        }           

        public void ReCompress()
        {
            // This is just a dumb algorithm that uses one tree entry per character and does not attempt compression.
            List<uint> newTreeEntries = new List<uint>();
            newTreeEntries.Add(0); // Dummy
            Dictionary<char, uint> entryForChar = new Dictionary<char, uint>();
            foreach (Id id in AllIds())
            {
                List<uint> newTreePointers = new List<uint>();
                foreach (char c in id.str)
                {
                    uint pointer;
                    if (!entryForChar.TryGetValue(c, out pointer))
                    {
                        pointer = (uint)newTreeEntries.Count;
                        newTreeEntries.Add(c);
                        entryForChar[c] = pointer;
                    }
                    newTreePointers.Add(pointer);
                }
                id.tree_pointers = newTreePointers.ToArray();
            }
            newTreeEntries[0] = (uint)newTreeEntries.Count;
            tree_entries = newTreeEntries.ToArray();
            tree_meta = new uint[] { 0xFFFFFFFF, 0x8 };
        }

    }

    public class Table
    {
        public uint first_id;
        public uint offset;

        public SubTableMeta[] sub_table_metas;
        public SubTableIds[] sub_table_ids;

        public long Read(Stream input, MemoryMappedViewAccessor accessor, Loc loc)
        {
            first_id = input.ReadValueU32();
            uint offset_length = input.ReadValueU32();

            long save = input.Position;

            offset = offset_length >> 4;
            uint length = offset_length & 15;

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

                sub_table_ids[i] = new SubTableIds();
                sub_table_ids[i].Read(ref block_first_id, sub_table_metas[i], input, accessor, loc);

                block_position += sub_table_metas[i].size;
            }

            input.Position = save;

            return block_position;
        }

        public void Write(Stream output)
        {
            output.WriteValueU32(first_id);
            output.WriteValueU32(offset << 4 | (uint)sub_table_metas.Length);
        }

        public List<Table> WriteData(Stream output, Loc loc)
        {
            List<SubTableMeta> newSubTableMetas = new List<SubTableMeta>();
            List<byte[]> allChunks = new List<byte[]>();
            uint id = first_id;
            for (int i = 0; i < sub_table_metas.Length; i++)
            {
                id += sub_table_metas[i].delta_from_prev_id;
                List<byte[]> chunks = sub_table_ids[i].Write(sub_table_metas[i], loc, id, newSubTableMetas);
                allChunks.AddRange(chunks);
                id += sub_table_metas[i].max_id + 1;
            }

            List<Table> tables = new List<Table>();
            Table pendingTable = new Table();
            pendingTable.offset = (uint)output.Position;
            pendingTable.first_id = first_id;
            id = first_id;
            List<SubTableMeta> pendingMetas = new List<SubTableMeta>();
            List<byte[]> pendingChunks = new List<byte[]>();
            for (int i = 0; i < newSubTableMetas.Count; i++)
            {
                id += newSubTableMetas[i].delta_from_prev_id + newSubTableMetas[i].max_id + 1;
                pendingMetas.Add(newSubTableMetas[i]);
                pendingChunks.Add(allChunks[i]);
                if (pendingMetas.Count == 15)
                {
                    pendingTable.sub_table_metas = pendingMetas.ToArray();
                    pendingMetas = new List<SubTableMeta>();
                    foreach (SubTableMeta meta in pendingTable.sub_table_metas)
                    {
                        meta.Write(output);
                    }
                    foreach (byte[] chunk in pendingChunks)
                    {
                        output.WriteBytes(chunk);
                    }
                    pendingChunks = new List<byte[]>();
                    tables.Add(pendingTable);
                    pendingTable = new Table();
                    pendingTable.offset = (uint)output.Position;
                    pendingTable.first_id = id;
                }
            }
            if (pendingMetas.Count > 0)
            {
                pendingTable.sub_table_metas = pendingMetas.ToArray();
                foreach (SubTableMeta meta in pendingTable.sub_table_metas)
                {
                    meta.Write(output);
                }
                foreach (byte[] chunk in pendingChunks)
                {
                    output.WriteBytes(chunk);
                }
                tables.Add(pendingTable);
            }
            return tables;
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

        public void Write(Stream output)
        {
            if (max_id <= 0x7F && size <= 0x7FF && delta_from_prev_id <= 0xFFF)
            {
                uint tmp = (max_id << 23) | (size << 12) | delta_from_prev_id;
                output.WriteValueU16((ushort)(tmp >> 16));
                output.WriteValueU16((ushort)(tmp & 0xFFFF));
            }
            else
            {
                uint tmp = 0x80000000 | (max_id << 16) | (delta_from_prev_id & 0xFFFF);
                if (delta_from_prev_id > 0xFFFF)
                {
                    tmp |= 0x40000000;
                }
                output.WriteValueU16((ushort)(tmp >> 16));
                output.WriteValueU16((ushort)(tmp & 0xFFFF));
                output.WriteValueU16((ushort)size);
                if (delta_from_prev_id > 0xFFFF)
                {
                    output.WriteValueU16((ushort)(delta_from_prev_id >> 16));
                }
            }
        }
    }

    public class SubTableIds
    {
        public List<List<Id>> ids;

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

        private List<Id> MakeFlatIdList()
        {
            List<Id> result = new List<Id>();
            foreach (List<Id> ids_in_block in ids)
            {
                foreach (Id id in ids_in_block)
                {
                    if (!id.is_pseudo)
                    {
                        result.Add(id);
                    }
                }
            }
            return result;
        }

        public List<byte[]> Write(SubTableMeta subTable, Loc loc, uint startId, List<SubTableMeta> subTables)
        {
            const int sizeTreshold = 65536 - 4096;
            List<byte[]> result = new List<byte[]>();
            List<Id> flatIds = MakeFlatIdList();
            SubTableMeta pendingSubTable = new SubTableMeta();
            pendingSubTable.delta_from_prev_id = subTable.delta_from_prev_id;
            int idx = 0;
            while (idx < flatIds.Count)
            {
                uint last_written_offset = 0xFFFFFFFF;
                List<ushort> block_64ids_offsets = new List<ushort>();
                MemoryStream tmp = new MemoryStream();
                MemoryStream tmpBits = new MemoryStream();
                BitWriter bitWriter = new BitWriter(tmpBits, 0);
                for (; idx < flatIds.Count; idx++)
                {
                    Id id = flatIds[idx];
                    uint offset = id.id - startId;
                    if (last_written_offset != 0xFFFFFFFF && (last_written_offset >> 6) != (offset >> 6))
                    {
                        // crossed a block treshold
                        uint rounder = ~last_written_offset & 0x3f;
                        Id.Skip(rounder, tmp);
                        last_written_offset += rounder;
                        while (((last_written_offset + 1) >> 6) != (offset >> 6))
                        {
                            block_64ids_offsets.Add(0);
                            last_written_offset += 64;
                        }
                        bitWriter.Close();
                        tmp.WriteBytes(tmpBits.ToArray());
                        tmpBits = new MemoryStream();
                        bitWriter = new BitWriter(tmpBits, 0);
                        block_64ids_offsets.Add((ushort)tmp.Position);
                    }
                    if (tmp.Position + tmpBits.Position >= sizeTreshold)
                    {
                        break;
                    }
                    Id.Skip(offset - last_written_offset - 1, tmp);
                    last_written_offset = offset - 1;
                    int before = bitWriter.Position;
                    foreach (uint ptr in id.tree_pointers)
                    {
                        loc.WriteTreePosition(bitWriter, ptr);
                    }
                    id.increment = (uint)(bitWriter.Position - before);
                    id.Write(tmp);
                    last_written_offset++;
                }
                if (bitWriter.Position == 0 && tmp.Position == block_64ids_offsets.Last())
                {
                    block_64ids_offsets.RemoveAt(block_64ids_offsets.Count - 1);
                }
                MemoryStream output = new MemoryStream();
                foreach (ushort offset in block_64ids_offsets)
                {
                    output.WriteValueU16((ushort)(offset == 0 ? 0 : offset + block_64ids_offsets.Count * 2));
                }
                output.WriteBytes(tmp.ToArray());
                bitWriter.Close();
                output.WriteBytes(tmpBits.ToArray());
                result.Add(output.ToArray());
                pendingSubTable.max_id = last_written_offset;
                pendingSubTable.size = (uint)output.Position;
                if (pendingSubTable.size > 0xFFFF)
                {
                    Console.WriteLine("Subtable got too big!");
                    Environment.Exit(1);
                }
                subTables.Add(pendingSubTable);
                if (idx < flatIds.Count)
                {
                    pendingSubTable = new SubTableMeta();
                    pendingSubTable.delta_from_prev_id = flatIds[idx].id - startId - last_written_offset - 1;
                    startId = flatIds[idx].id;
                }
            }
            return result;
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

        public Id(bool is_pseudo, uint increment)
        {
            this.is_pseudo = is_pseudo;
            this.increment = increment;
        }

        public static void Skip(uint count, Stream output)
        {
            while (count > 0)
            {
                uint page = Math.Min(15, count);
                new Id(true, page).Write(output);
                count -= page;
            }
        }

        public override string ToString()
        {
            return "id=" + id + " increment=" + (is_pseudo ? "*" : "") + increment + " tree_pointers=" + tree_pointers + " str=" + str;
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

        public void Write(Stream output)
        {
            if (is_pseudo)
            {
                output.WriteValueU8((byte)(increment + 0xF0));
            }
            else
            {
                uint adjusted_size = increment == 0 ? 0 : (increment - 4) / 2;
                if (adjusted_size >= 0x114dc)
                {
                    Console.WriteLine("Bit string too big: " + adjusted_size);
                    Environment.Exit(2);
                }
                else if (adjusted_size >= 0x14dc)
                {
                    output.WriteValueU8(0xF0);
                    uint tmp = adjusted_size - 5340;
                    output.WriteValueU8((byte)(tmp >> 8));
                    output.WriteValueU8((byte)(tmp & 0xFF));
                }
                else if (adjusted_size >= 0xdc)
                {
                    uint tmp = adjusted_size + 56100;
                    output.WriteValueU8((byte)(tmp >> 8));
                    output.WriteValueU8((byte)(tmp & 0xFF));
                }
                else
                {
                    output.WriteValueU8((byte)adjusted_size);
                }
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

    public class BitWriter
    {
        readonly Stream stream;
        ulong pendingBits;
        int pendingBitCount;

        public BitWriter(Stream stream, int position)
        {
            this.stream = stream;
            this.Position = position;
        }

        public int Position { get; private set; }

        public void Write(uint bits, int bitCount)
        {
            if (bitCount < 1 || bitCount > 32)
            {
                throw new ArgumentException("Bit count out of range: " + bitCount);
            }
            pendingBits |= (ulong)bits << (64 - bitCount - pendingBitCount);
            pendingBitCount += bitCount;
            while (pendingBitCount >= 8)
            {
                stream.WriteValueU8((byte)(pendingBits >> 56));
                pendingBits <<= 8;
                pendingBitCount -= 8;
            }
            Position += bitCount;
        }

        public void Close()
        {
            if (pendingBitCount > 0)
            {
                stream.WriteValueU8((byte)(pendingBits >> 56));
            }
        }
    }
}
