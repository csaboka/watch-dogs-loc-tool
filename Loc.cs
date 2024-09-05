using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
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

        public void DecodeStrings()
        {
            foreach (Id id in AllIds())
            {
                id.str = DecodeString(id.tree_pointers);
            }
        }

        public void Export(String filename)
        {
            using (StreamWriter text = new StreamWriter(filename + ".txt", false, System.Text.Encoding.Unicode))
            {
                foreach (Id id in AllIds())
                {
                    text.WriteLine(id.id + "=" + id.str.Replace("\r", "[CR]").Replace("\n", "[LF]"));
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
            HashSet<uint> unseenKeys = new HashSet<uint>(newStrings.Keys);
            foreach (Id id in AllIds())
            {
                string newValue;
                if (newStrings.TryGetValue(id.id, out newValue))
                {
                    id.str = newValue;
                    unseenKeys.Remove(id.id);
                } else
                {
                    Console.WriteLine("WARNING: ID {0} missing from the text file! Will keep current value from the LOC file", id.id);
                }
            }
            if (unseenKeys.Count > 0)
            {
                Console.WriteLine("WARNING: IDs {0} were present in the text file but not in the LOC file. They were ignored. Please check the text file for typos.", string.Join(", ", unseenKeys));
            }
        }

        public void ReCompress(bool highEffort)
        {
            Stopwatch grandTotal = new Stopwatch();
            grandTotal.Start();
            MakeStarterEntries();
            if (highEffort)
            {
                // Until I figure out how to  rebuild tree_meta in the generic case, it isn't safe to use tree indices above
                // what the original file could encode.
                AddCompressionEntries((int)(((tree_meta[11] >> 8) - 1) & 0xFF00));
                ShuffleEntriesByFrequency();
            }

            tree_entries[0] = (uint)tree_entries.Length;
            if (tree_entries.Length < 0xFF)
            {
                // Looks like the game expects 12 tree_meta entries, so build a fake table where the first entry is always chosen.
                // To be on the safe side, also re-use the original bit lengths.
                tree_meta = new uint[] {
                    0xFF000000, 0x8,
                    0xFF010000, 0xa,
                    0xFF020000, 0xc,
                    0xFF030000, 0xe,
                    0xFF040000, 0x10,
                    0xFFFFFFFF, 0x18
                };
            }
            if (highEffort)
            {
                Console.WriteLine("Total time spent compressing: " + grandTotal.Elapsed);
            }
        }

        void MakeStarterEntries()
        {
            List<uint> newTreeEntries = new List<uint>();
            newTreeEntries.Add(0); // Dummy
            Dictionary<char, uint> entryForChar = new Dictionary<char, uint>();
            foreach (Id id in AllIds())
            {
                List<uint> newTreePointers = new List<uint>(id.str.Length);
                foreach (char c in id.str)
                {
                    uint pointer;
                    if (!entryForChar.TryGetValue(c, out pointer))
                    {
                        pointer = (ushort)newTreeEntries.Count;
                        newTreeEntries.Add(c);
                        entryForChar[c] = pointer;
                    }
                    newTreePointers.Add(pointer);
                }
                id.tree_pointers = newTreePointers.ToArray();
            }
            tree_entries = newTreeEntries.ToArray();
        }

        static long BloomFor(ushort value)
        {
            return 1L << (value & 63);
        }

        void AddCompressionEntries(int maxEntries)
        {
            // The general idea is that to optimize the file, we pick the pair of symbols that appear most frequently, introduce
            // a new symbol for them, then replace all occurrences of the pair with the new symbol. Then we repeat this process,
            // generating a new symbol and eliminating some tree pointers every time, until we run out of symbols we can allocate.

            // This process requires iterating through the entries an ungodly amount of times, so we do our best to optimize things:
            // - Pointers are converted from uints to ushorts to improve cache locality.
            // - We build a simple 64-bit Bloom filter for each entry. This lets us quickly skip entries that don't contain both
            //   symbols of the current pair. This starts out somewhat poorly (around 50% passing through, around 50% false positives),
            //   but improves drastically after roughly a thousand iterations.
            // - We periodically trim the pointer lists to bring the remaining items closer together in memory and improve cache locality.
            List<uint> newTreeEntries = new List<uint>(tree_entries);
            List<List<ushort>> allPointersList = new List<List<ushort>>();
            PairFrequencies pairFrequencies = new PairFrequencies();
            int total_len = 0;
            List<long> bloomList = new List<long>();
            foreach (Id id in AllIds())
            {
                total_len += id.tree_pointers.Length;
                List<ushort> shortList = new List<ushort>(id.tree_pointers.Length);
                shortList.AddRange(id.tree_pointers.Select(ptr => (ushort)ptr));
                allPointersList.Add(shortList);
                for (int i = 0; i < shortList.Count - 1; i++)
                {
                    pairFrequencies.AddPair(shortList[i], shortList[i + 1]);
                }
                long bloomThis = 0;
                foreach (ushort pointer in shortList)
                {
                    bloomThis |= BloomFor(pointer);
                }
                bloomList.Add(bloomThis);
            }
            int starting_total_len = total_len;

            List<ushort>[] allPointersArray = allPointersList.ToArray();
            long[] bloomArray = bloomList.ToArray();
            // Iterating through the frequency table is pretty expensive, so in addition to picking a winner,
            // collect all the ties as well. Then consume all the ties before going through the frequency table again.
            // We need to revalidate the pending items because their frequency may have been decreased by the previous
            // replacements, but even with that step, this approach is much faster than a naive one.
            Queue<PairFrequency> tiedForBest = new Queue<PairFrequency>();
            while (newTreeEntries.Count < maxEntries)
            {
                PairFrequency bestPair;
                while (tiedForBest.Count > 0)
                {
                    PairFrequency top = tiedForBest.Peek();
                    if (top.count == pairFrequencies.GetPairFrequency(top.first, top.second))
                    {
                        break;
                    }
                    else
                    {
                        // stale entry, get rid of it
                        tiedForBest.Dequeue();
                    }
                }
                if (tiedForBest.Count > 0)
                {
                    bestPair = tiedForBest.Dequeue();
                }
                else
                {
                    bestPair = new PairFrequency();
                    foreach (PairFrequency pair in pairFrequencies)
                    {
                        if (pair.count > bestPair.count)
                        {
                            bestPair = pair;
                            tiedForBest.Clear();
                        }
                        else if (pair.count == bestPair.count)
                        {
                            tiedForBest.Enqueue(pair);
                        }
                    }
                }
                if ((newTreeEntries.Count & 0xFF) == 0)
                {
                    if ((newTreeEntries.Count & 2047) == 0)
                    {
                        foreach (List<ushort> pointers in allPointersArray)
                        {
                            pointers.TrimExcess();
                        }
                    }
                    Console.Write("Progress: {0:P}, estimated size factor: {1:P}\r", (double)newTreeEntries.Count / maxEntries, (double)total_len / starting_total_len);
                }
                ushort pointer = (ushort)newTreeEntries.Count;
                newTreeEntries.Add((uint)bestPair.first << 16 | bestPair.second);
                long filter = BloomFor(bestPair.first) | BloomFor(bestPair.second);
                // replace all pairs with this new entry while updating the bookkeeping
                for (int bloomIdx = 0; bloomIdx < bloomArray.Length; bloomIdx++)
                {
                    if ((bloomArray[bloomIdx] & filter) != filter)
                    {
                        continue;
                    }
                    List<ushort> pointers = allPointersArray[bloomIdx];
                    int startIndex = 0;
                    bool hadChange = false;
                    while (startIndex < pointers.Count)
                    {
                        int idx = pointers.IndexOf(bestPair.first, startIndex);
                        if (idx < 0 || idx == pointers.Count - 1)
                        {
                            break;
                        }
                        if (pointers[idx + 1] == bestPair.second)
                        {
                            if (idx > 0)
                            {
                                pairFrequencies.RemovePair(pointers[idx - 1], pointers[idx]);
                                pairFrequencies.AddPair(pointers[idx - 1], pointer);
                            }
                            pairFrequencies.RemovePair(bestPair.first, bestPair.second);
                            if (idx < pointers.Count - 2)
                            {
                                pairFrequencies.RemovePair(pointers[idx + 1], pointers[idx + 2]);
                                pairFrequencies.AddPair(pointer, pointers[idx + 2]);
                            }
                            pointers.RemoveAt(idx + 1);
                            total_len--;
                            pointers[idx] = pointer;
                            hadChange = true;
                        }
                        startIndex = idx + 1;
                    }
                    if (hadChange)
                    {
                        bloomArray[bloomIdx] = 0;
                        foreach (ushort ptr in pointers)
                        {
                            bloomArray[bloomIdx] |= BloomFor(ptr);
                        }
                    }
                }
            }

            int repackIndex = 0;
            foreach (Id id in AllIds())
            {
                id.tree_pointers = allPointersArray[repackIndex++].Select(ptr => (uint)ptr).ToArray();
            }
            tree_entries = newTreeEntries.ToArray();

            Console.WriteLine("Progress: {0:P}, estimated size factor: {1:P}", 1, (double)total_len / starting_total_len);
        }

        void ShuffleEntriesByFrequency()
        {
            // Lower indices are cheaper to encode, so reorder the tree to make the most popular symbols appear first.
            // This will also neatly move the symbols that aren't used directly (i.e. only referenced by other tree entries)
            // to the very end. Since those entries never need to be encoded, any index is fine for them.
            uint[] frequencies = new uint[tree_entries.Length];
            foreach (Id id in AllIds())
            {
                foreach (uint ptr in id.tree_pointers)
                {
                    frequencies[ptr]++;
                }
            }
            // Index 0 is special, make sure it stays the first by giving it a dummy frequency.
            frequencies[0] = uint.MaxValue;
            uint[] translationTable = new uint[tree_entries.Length];
            uint newIdx = 0;
            foreach (int oldIdx in Enumerable.Range(0, tree_entries.Length).OrderByDescending(i => frequencies[i]))
            {
                translationTable[oldIdx] = newIdx++;
            }

            foreach (Id id in AllIds())
            {
                id.tree_pointers = id.tree_pointers.Select(w => translationTable[w]).ToArray();
            }

            uint[] newTreeEntries = new uint[tree_entries.Length];
            for (int i=1; i< tree_entries.Length; i++)
            {
                uint entry = tree_entries[i];
                if (entry > 0xFFFF)
                {
                    uint first = translationTable[entry >> 16];
                    uint second = translationTable[entry & 0xFFFF];
                    entry = (first << 16) | second;
                }
                newTreeEntries[translationTable[i]] = entry;
            }
            tree_entries = newTreeEntries;
        }

    }

    struct PairFrequency
    {
        public ushort first;
        public ushort second;
        public uint count;
    }

    class PairFrequencies : IEnumerable<PairFrequency>
    {
        private Dictionary<uint, uint> freq = new Dictionary<uint, uint>();

        uint MakeKey(ushort first, ushort second)
        {
            return ((uint)first << 16) | second;
        }

        public void AddPair(ushort first, ushort second)
        {
            uint key = MakeKey(first, second);
            uint count;
            freq.TryGetValue(key, out count);
            freq[key] = count + 1;
        }

        public void RemovePair(ushort first, ushort second)
        {
            uint key = MakeKey(first, second);
            uint newCount = --freq[key];
            if (newCount == 0)
            {
                freq.Remove(key);
            }
        }

        public uint GetPairFrequency(ushort first, ushort second)
        {
            freq.TryGetValue(MakeKey(first, second), out uint result);
            return result;
        }

        public IEnumerator<PairFrequency> GetEnumerator()
        {
            foreach (var pair in freq)
            {
                yield return new PairFrequency
                {
                    first = (ushort)(pair.Key >> 16),
                    second = (ushort)(pair.Key & 0xFFFF),
                    count = pair.Value
                };
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            foreach (PairFrequency pf in this)
            {
                yield return pf;
            }
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
                bool largeDelta = delta_from_prev_id > 0xFFFF;
                uint tmp = (largeDelta ? 0xc0000000 : 0x80000000) | (max_id << 16) | (delta_from_prev_id & 0xFFFF);
                output.WriteValueU16((ushort)(tmp >> 16));
                output.WriteValueU16((ushort)(tmp & 0xFFFF));
                output.WriteValueU16((ushort)size);
                if (largeDelta)
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
                if (j >= 64)
                {
                    long block_64ids_offset = block_64ids_offsets[(j >> 6) - 1];
                    if (block_64ids_offset == 0)
                    {
                        ids.Add(ids_in_block);
                        continue;
                    }
                    input.Position = subtable_ids_begin + block_64ids_offset;
                }

               
                uint k = 0;
                while (k < Math.Min(id_count - j, 64))
                {
                    Id new_id = new Id(id_begin + j + k);
                    new_id.Read(ref k, input);
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
                BitWriter bitWriter = new BitWriter();
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
                        tmp.WriteBytes(bitWriter.Finish());
                        bitWriter = new BitWriter();
                        block_64ids_offsets.Add((ushort)tmp.Position);
                    }
                    if (tmp.Position + (bitWriter.Position / 8) >= sizeTreshold)
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
                output.WriteBytes(bitWriter.Finish());
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

        public void Read(ref uint k, Stream input)
        {
            uint current_size = input.ReadValueU8();
            if (current_size > 0xF0)
            {
                increment = current_size - 240;
                k += increment;
                if (k > 64)
                {
                    Console.WriteLine("Extras! last id=" + id + " skip=" + increment + " extra=" + (k - 64));
                    Environment.Exit(1);
                }
                is_pseudo = true;
            }
            else
            {
                if (current_size == 0xF0)
                {
                    current_size = input.ReadValueU8();
                    current_size = (current_size << 8) + input.ReadValueU8() + 5340;
                }
                else if (current_size >= 0xDC)
                {
                    current_size = (current_size << 8) + input.ReadValueU8() - 56100;
                }
                if (current_size != 0)
                {
                    current_size = 2 * current_size + 4;
                }
                increment = current_size;
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
        List<byte> buffer = new List<byte>();
        ulong pendingBits;
        int pendingBitCount;

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
                buffer.Add((byte)(pendingBits >> 56));
                pendingBits <<= 8;
                pendingBitCount -= 8;
            }
            Position += bitCount;
        }

        public byte[] Finish()
        {
            if (pendingBitCount > 0)
            {
                buffer.Add((byte)(pendingBits >> 56));
            }
            return buffer.ToArray();
        }
    }
}
