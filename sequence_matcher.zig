const std = @import("std");
const mem = std.mem;
const assert = std.debug.assert;
const ArrayList = std.ArrayList;
const AutoHashMap = std.AutoHashMap;

pub const LongestMatch = struct {
    a_idx: usize,
    b_idx: usize,
    size: usize,
};

pub const SequenceMatcher = struct {
    seq1: []const u8,
    seq2: []const u8,
    allocator: *mem.Allocator,

    b2j: AutoHashMap(u8, ArrayList(usize)),
    matching_blocks: []i32, 
    opcodes: []i32,
    full_b_count: []i32,

    pub fn init(allocator: *mem.Allocator, a: []const u8, b: []const u8) !SequenceMatcher {
        var b2j = AutoHashMap(u8, ArrayList(usize)).init(allocator);
        var sm = SequenceMatcher{
            .allocator = allocator,
            .seq1 = undefined,
            .seq2 = undefined,
            .b2j = b2j,
            .matching_blocks = undefined,
            .opcodes = undefined,
            .full_b_count = undefined,
        };
        try sm.setSeqs(a, b);
        return sm;
    }

    pub fn deinit(self: *SequenceMatcher) void {
        var it = self.b2j.iterator();
        while (it.next()) |indices| {
            indices.value.deinit();
        }
        self.b2j.deinit();
    }
    
    pub fn setSeqs(self: *SequenceMatcher, a: []const u8, b: []const u8) !void {
        self.setSeq1(a);
        try self.setSeq2(b);
    }

    pub fn setSeq1(self: *SequenceMatcher, a: []const u8) void {
        if (mem.eql(u8, a, self.seq1)) {
            return;
        }

        self.seq1 = a;
        self.matching_blocks = undefined;
        self.opcodes = undefined;
    }

    pub fn setSeq2(self: *SequenceMatcher, b: []const u8) !void {
        if (mem.eql(u8, b, self.seq2)) {
            return;
        }
        self.seq2 = b;
        self.matching_blocks = undefined;
        self.opcodes = undefined;
        self.full_b_count = undefined;

        try self.chainB();
    }

    fn chainB(self: *SequenceMatcher) !void {
        var b = self.seq2;
        self.b2j.clear();

        var n_test = @intToFloat(f32, b.len) / 100 + 1;

        for (b) |byte, idx| {
            if (self.b2j.getValue(byte)) |*byte_idxs| {
                var byte_idx_size = @intToFloat(f32, byte_idxs.items.len);
                if (b.len >= 200 and byte_idx_size > n_test) continue;
                try byte_idxs.append(idx);
                _ = try self.b2j.put(byte, byte_idxs.*);
            } else {
                var byte_idxs = ArrayList(usize).init(self.allocator);
                try byte_idxs.append(idx);
                _ = try self.b2j.put(byte, byte_idxs);
            }
        }
    }

    fn findLongestMatch(self: *SequenceMatcher, a_lo: usize, a_hi: usize, b_lo: usize, b_hi: usize) !LongestMatch {
        var a = self.seq1;
        var b = self.seq2;
        assert(a_hi <= a.len and b_hi <= b.len);

        var b2j = self.b2j;
        var best_i = a_lo;
        var best_j = b_lo;
        var best_size: usize = 0;
        var j2_len = AutoHashMap(usize, usize).init(self.allocator);
        var nothing = [0]usize{};

        var i = a_lo;
        while (i < a_hi) : (i += 1) {
            var indices: []usize = undefined;
            if (b2j.getValue(a[i])) |value| {
                indices = value.items;
            } else {
                indices = &nothing;
            }

            var next_j2_len: AutoHashMap(usize, usize) = undefined;
            var updated_j2_len = false;
            for (indices) |j| {
                if (j < b_lo) continue;
                if (j >= b_hi) break;
                var k = (if (j == 0) 0 else (j2_len.getValue(j - 1) orelse 0)) + 1;

                if (!updated_j2_len) {
                    next_j2_len = AutoHashMap(usize, usize).init(self.allocator);
                    updated_j2_len = true;
                }
                _ = try next_j2_len.put(j, k);
                if (k > best_size) {
                    best_i = i + 1 - k;
                    best_j = j + 1 - k;
                    best_size = k;
                }
            }
            if (updated_j2_len) {
                j2_len.deinit();
                j2_len = next_j2_len;
            }
        }
        j2_len.deinit();

        while (best_i > a_lo and best_j > b_lo and a[best_i - 1] == b[best_j - 1]) {
            best_i -= 1;
            best_j -= 1;
            best_size += 1;
        }

        while (best_i + best_size < a_hi and best_j + best_size < b_hi and a[best_i + best_size] == b[best_j + best_size]) {
            best_size += 1;
        }

        return LongestMatch{
            .a_idx = best_i,
            .b_idx = best_j,
            .size = best_size,
        };
    }
};

test "sequence matcher - find longest match" {
    const allocator = std.testing.allocator;
    // var sm = try SequenceMatcher.init(allocator, "abcd", "bcde");
    var sm = try SequenceMatcher.init(allocator, " abcd", "abcd abcd");
    defer sm.deinit();
    var lm = try sm.findLongestMatch(0, 5, 0, 9);
    assert(lm.a_idx == 0);
    assert(lm.b_idx == 4);
    assert(lm.size == 5);
}

//     def get_matching_blocks(self):
//         """Return list of triples describing matching subsequences.
//         Each triple is of the form (i, j, n), and means that
//         a[i:i+n] == b[j:j+n].  The triples are monotonically increasing in
//         i and in j.  New in Python 2.5, it's also guaranteed that if
//         (i, j, n) and (i', j', n') are adjacent triples in the list, and
//         the second is not the last triple in the list, then i+n != i' or
//         j+n != j'.  IOW, adjacent triples never describe adjacent equal
//         blocks.
//         The last triple is a dummy, (len(a), len(b), 0), and is the only
//         triple with n==0.
//         >>> s = SequenceMatcher(None, "abxcd", "abcd")
//         >>> list(s.get_matching_blocks())
//         [Match(a=0, b=0, size=2), Match(a=3, b=2, size=2), Match(a=5, b=4, size=0)]
//         """

//         if self.matching_blocks is not None:
//             return self.matching_blocks
//         la, lb = len(self.a), len(self.b)

//         # This is most naturally expressed as a recursive algorithm, but
//         # at least one user bumped into extreme use cases that exceeded
//         # the recursion limit on their box.  So, now we maintain a list
//         # ('queue`) of blocks we still need to look at, and append partial
//         # results to `matching_blocks` in a loop; the matches are sorted
//         # at the end.
//         queue = [(0, la, 0, lb)]
//         matching_blocks = []
//         while queue:
//             alo, ahi, blo, bhi = queue.pop()
//             i, j, k = x = self.find_longest_match(alo, ahi, blo, bhi)
//             # a[alo:i] vs b[blo:j] unknown
//             # a[i:i+k] same as b[j:j+k]
//             # a[i+k:ahi] vs b[j+k:bhi] unknown
//             if k:   # if k is 0, there was no matching block
//                 matching_blocks.append(x)
//                 if alo < i and blo < j:
//                     queue.append((alo, i, blo, j))
//                 if i+k < ahi and j+k < bhi:
//                     queue.append((i+k, ahi, j+k, bhi))
//         matching_blocks.sort()

//         # It's possible that we have adjacent equal blocks in the
//         # matching_blocks list now.  Starting with 2.5, this code was added
//         # to collapse them.
//         i1 = j1 = k1 = 0
//         non_adjacent = []
//         for i2, j2, k2 in matching_blocks:
//             # Is this block adjacent to i1, j1, k1?
//             if i1 + k1 == i2 and j1 + k1 == j2:
//                 # Yes, so collapse them -- this just increases the length of
//                 # the first block by the length of the second, and the first
//                 # block so lengthened remains the block to compare against.
//                 k1 += k2
//             else:
//                 # Not adjacent.  Remember the first block (k1==0 means it's
//                 # the dummy we started with), and make the second block the
//                 # new block to compare against.
//                 if k1:
//                     non_adjacent.append((i1, j1, k1))
//                 i1, j1, k1 = i2, j2, k2
//         if k1:
//             non_adjacent.append((i1, j1, k1))

//         non_adjacent.append( (la, lb, 0) )
//         self.matching_blocks = list(map(Match._make, non_adjacent))
//         return self.matching_blocks

//     def get_opcodes(self):
//         """Return list of 5-tuples describing how to turn a into b.
//         Each tuple is of the form (tag, i1, i2, j1, j2).  The first tuple
//         has i1 == j1 == 0, and remaining tuples have i1 == the i2 from the
//         tuple preceding it, and likewise for j1 == the previous j2.
//         The tags are strings, with these meanings:
//         'replace':  a[i1:i2] should be replaced by b[j1:j2]
//         'delete':   a[i1:i2] should be deleted.
//                     Note that j1==j2 in this case.
//         'insert':   b[j1:j2] should be inserted at a[i1:i1].
//                     Note that i1==i2 in this case.
//         'equal':    a[i1:i2] == b[j1:j2]
//         >>> a = "qabxcd"
//         >>> b = "abycdf"
//         >>> s = SequenceMatcher(None, a, b)
//         >>> for tag, i1, i2, j1, j2 in s.get_opcodes():
//         ...    print(("%7s a[%d:%d] (%s) b[%d:%d] (%s)" %
//         ...           (tag, i1, i2, a[i1:i2], j1, j2, b[j1:j2])))
//          delete a[0:1] (q) b[0:0] ()
//           equal a[1:3] (ab) b[0:2] (ab)
//         replace a[3:4] (x) b[2:3] (y)
//           equal a[4:6] (cd) b[3:5] (cd)
//          insert a[6:6] () b[5:6] (f)
//         """

//         if self.opcodes is not None:
//             return self.opcodes
//         i = j = 0
//         self.opcodes = answer = []
//         for ai, bj, size in self.get_matching_blocks():
//             # invariant:  we've pumped out correct diffs to change
//             # a[:i] into b[:j], and the next matching block is
//             # a[ai:ai+size] == b[bj:bj+size].  So we need to pump
//             # out a diff to change a[i:ai] into b[j:bj], pump out
//             # the matching block, and move (i,j) beyond the match
//             tag = ''
//             if i < ai and j < bj:
//                 tag = 'replace'
//             elif i < ai:
//                 tag = 'delete'
//             elif j < bj:
//                 tag = 'insert'
//             if tag:
//                 answer.append( (tag, i, ai, j, bj) )
//             i, j = ai+size, bj+size
//             # the list of matching blocks is terminated by a
//             # sentinel with size 0
//             if size:
//                 answer.append( ('equal', ai, i, bj, j) )
//         return answer

//     def get_grouped_opcodes(self, n=3):
//         """ Isolate change clusters by eliminating ranges with no changes.
//         Return a generator of groups with up to n lines of context.
//         Each group is in the same format as returned by get_opcodes().
//         >>> from pprint import pprint
//         >>> a = list(map(str, range(1,40)))
//         >>> b = a[:]
//         >>> b[8:8] = ['i']     # Make an insertion
//         >>> b[20] += 'x'       # Make a replacement
//         >>> b[23:28] = []      # Make a deletion
//         >>> b[30] += 'y'       # Make another replacement
//         >>> pprint(list(SequenceMatcher(None,a,b).get_grouped_opcodes()))
//         [[('equal', 5, 8, 5, 8), ('insert', 8, 8, 8, 9), ('equal', 8, 11, 9, 12)],
//          [('equal', 16, 19, 17, 20),
//           ('replace', 19, 20, 20, 21),
//           ('equal', 20, 22, 21, 23),
//           ('delete', 22, 27, 23, 23),
//           ('equal', 27, 30, 23, 26)],
//          [('equal', 31, 34, 27, 30),
//           ('replace', 34, 35, 30, 31),
//           ('equal', 35, 38, 31, 34)]]
//         """

//         codes = self.get_opcodes()
//         if not codes:
//             codes = [("equal", 0, 1, 0, 1)]
//         # Fixup leading and trailing groups if they show no changes.
//         if codes[0][0] == 'equal':
//             tag, i1, i2, j1, j2 = codes[0]
//             codes[0] = tag, max(i1, i2-n), i2, max(j1, j2-n), j2
//         if codes[-1][0] == 'equal':
//             tag, i1, i2, j1, j2 = codes[-1]
//             codes[-1] = tag, i1, min(i2, i1+n), j1, min(j2, j1+n)

//         nn = n + n
//         group = []
//         for tag, i1, i2, j1, j2 in codes:
//             # End the current group and start a new one whenever
//             # there is a large range with no changes.
//             if tag == 'equal' and i2-i1 > nn:
//                 group.append((tag, i1, min(i2, i1+n), j1, min(j2, j1+n)))
//                 yield group
//                 group = []
//                 i1, j1 = max(i1, i2-n), max(j1, j2-n)
//             group.append((tag, i1, i2, j1 ,j2))
//         if group and not (len(group)==1 and group[0][0] == 'equal'):
//             yield group

//     def ratio(self):
//         """Return a measure of the sequences' similarity (float in [0,1]).
//         Where T is the total number of elements in both sequences, and
//         M is the number of matches, this is 2.0*M / T.
//         Note that this is 1 if the sequences are identical, and 0 if
//         they have nothing in common.
//         .ratio() is expensive to compute if you haven't already computed
//         .get_matching_blocks() or .get_opcodes(), in which case you may
//         want to try .quick_ratio() or .real_quick_ratio() first to get an
//         upper bound.
//         >>> s = SequenceMatcher(None, "abcd", "bcde")
//         >>> s.ratio()
//         0.75
//         >>> s.quick_ratio()
//         0.75
//         >>> s.real_quick_ratio()
//         1.0
//         """

//         matches = sum(triple[-1] for triple in self.get_matching_blocks())
//         return _calculate_ratio(matches, len(self.a) + len(self.b))

//     def quick_ratio(self):
//         """Return an upper bound on ratio() relatively quickly.
//         This isn't defined beyond that it is an upper bound on .ratio(), and
//         is faster to compute.
//         """

//         # viewing a and b as multisets, set matches to the cardinality
//         # of their intersection; this counts the number of matches
//         # without regard to order, so is clearly an upper bound
//         if self.fullbcount is None:
//             self.fullbcount = fullbcount = {}
//             for elt in self.b:
//                 fullbcount[elt] = fullbcount.get(elt, 0) + 1
//         fullbcount = self.fullbcount
//         # avail[x] is the number of times x appears in 'b' less the
//         # number of times we've seen it in 'a' so far ... kinda
//         avail = {}
//         availhas, matches = avail.__contains__, 0
//         for elt in self.a:
//             if availhas(elt):
//                 numb = avail[elt]
//             else:
//                 numb = fullbcount.get(elt, 0)
//             avail[elt] = numb - 1
//             if numb > 0:
//                 matches = matches + 1
//         return _calculate_ratio(matches, len(self.a) + len(self.b))

//     def real_quick_ratio(self):
//         """Return an upper bound on ratio() very quickly.
//         This isn't defined beyond that it is an upper bound on .ratio(), and
//         is faster to compute than either .ratio() or .quick_ratio().
//         """

//         la, lb = len(self.a), len(self.b)
//         # can't have more matches than the number of elements in the
//         # shorter sequence
//         return _calculate_ratio(min(la, lb), la + lb)