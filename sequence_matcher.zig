const std = @import("std");
const mem = std.mem;
const unicode = std.unicode;
const ArrayList = std.ArrayList;
const AutoHashMap = std.AutoHashMap;

pub const SequenceMatcher = struct {
    seq1: []const u8,
    seq2: []const u8,
    allocator: *mem.Allocator,

    b2j: AutoHashMap(u21, ArrayList(usize)),
    matching_blocks: []i32, 
    opcodes: []i32,
    full_b_count: []i32,

    pub fn init(allocator: *mem.Allocator, a: []const u8, b: []const u8) !SequenceMatcher {
        var b2j = AutoHashMap(u21, ArrayList(usize)).init(allocator);
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
        var b2f = self.b2j;

        var n_test = @intToFloat(f32, b.len) / 100 + 1;

        var idx: usize = 0;
        var utf8 = (try unicode.Utf8View.init(b)).iterator();
        while (utf8.nextCodepoint()) |codepoint| {
            if (self.b2j.getValue(codepoint)) |*codepoint_idxs| {
                var codepoint_size = @intToFloat(f32, codepoint_idxs.items.len);
                if (b.len >= 200 and codepoint_size > n_test) {
                    continue;
                }
                try codepoint_idxs.append(idx);
            } else {
                var codepoint_idxs = ArrayList(usize).init(self.allocator);
                try codepoint_idxs.append(codepoint);
                _ = try self.b2j.put(codepoint, codepoint_idxs);
            }
        }

    }
};

test "sequence matcher" {
    const allocator = std.testing.allocator;
    var sm = try SequenceMatcher.init(allocator, "abcd", "bcde");
    defer sm.deinit();
}

//     def find_longest_match(self, alo, ahi, blo, bhi):
//         """Find longest matching block in a[alo:ahi] and b[blo:bhi].
//         If isjunk is not defined:
//         Return (i,j,k) such that a[i:i+k] is equal to b[j:j+k], where
//             alo <= i <= i+k <= ahi
//             blo <= j <= j+k <= bhi
//         and for all (i',j',k') meeting those conditions,
//             k >= k'
//             i <= i'
//             and if i == i', j <= j'
//         In other words, of all maximal matching blocks, return one that
//         starts earliest in a, and of all those maximal matching blocks that
//         start earliest in a, return the one that starts earliest in b.
//         >>> s = SequenceMatcher(None, " abcd", "abcd abcd")
//         >>> s.find_longest_match(0, 5, 0, 9)
//         Match(a=0, b=4, size=5)
//         If isjunk is defined, first the longest matching block is
//         determined as above, but with the additional restriction that no
//         junk element appears in the block.  Then that block is extended as
//         far as possible by matching (only) junk elements on both sides.  So
//         the resulting block never matches on junk except as identical junk
//         happens to be adjacent to an "interesting" match.
//         Here's the same example as before, but considering blanks to be
//         junk.  That prevents " abcd" from matching the " abcd" at the tail
//         end of the second sequence directly.  Instead only the "abcd" can
//         match, and matches the leftmost "abcd" in the second sequence:
//         >>> s = SequenceMatcher(lambda x: x==" ", " abcd", "abcd abcd")
//         >>> s.find_longest_match(0, 5, 0, 9)
//         Match(a=1, b=0, size=4)
//         If no blocks match, return (alo, blo, 0).
//         >>> s = SequenceMatcher(None, "ab", "c")
//         >>> s.find_longest_match(0, 2, 0, 1)
//         Match(a=0, b=0, size=0)
//         """

//         # CAUTION:  stripping common prefix or suffix would be incorrect.
//         # E.g.,
//         #    ab
//         #    acab
//         # Longest matching block is "ab", but if common prefix is
//         # stripped, it's "a" (tied with "b").  UNIX(tm) diff does so
//         # strip, so ends up claiming that ab is changed to acab by
//         # inserting "ca" in the middle.  That's minimal but unintuitive:
//         # "it's obvious" that someone inserted "ac" at the front.
//         # Windiff ends up at the same place as diff, but by pairing up
//         # the unique 'b's and then matching the first two 'a's.

//         a, b, b2j, isbjunk = self.a, self.b, self.b2j, self.bjunk.__contains__
//         besti, bestj, bestsize = alo, blo, 0
//         # find longest junk-free match
//         # during an iteration of the loop, j2len[j] = length of longest
//         # junk-free match ending with a[i-1] and b[j]
//         j2len = {}
//         nothing = []
//         for i in range(alo, ahi):
//             # look at all instances of a[i] in b; note that because
//             # b2j has no junk keys, the loop is skipped if a[i] is junk
//             j2lenget = j2len.get
//             newj2len = {}
//             for j in b2j.get(a[i], nothing):
//                 # a[i] matches b[j]
//                 if j < blo:
//                     continue
//                 if j >= bhi:
//                     break
//                 k = newj2len[j] = j2lenget(j-1, 0) + 1
//                 if k > bestsize:
//                     besti, bestj, bestsize = i-k+1, j-k+1, k
//             j2len = newj2len

//         # Extend the best by non-junk elements on each end.  In particular,
//         # "popular" non-junk elements aren't in b2j, which greatly speeds
//         # the inner loop above, but also means "the best" match so far
//         # doesn't contain any junk *or* popular non-junk elements.
//         while besti > alo and bestj > blo and \
//               not isbjunk(b[bestj-1]) and \
//               a[besti-1] == b[bestj-1]:
//             besti, bestj, bestsize = besti-1, bestj-1, bestsize+1
//         while besti+bestsize < ahi and bestj+bestsize < bhi and \
//               not isbjunk(b[bestj+bestsize]) and \
//               a[besti+bestsize] == b[bestj+bestsize]:
//             bestsize += 1

//         # Now that we have a wholly interesting match (albeit possibly
//         # empty!), we may as well suck up the matching junk on each
//         # side of it too.  Can't think of a good reason not to, and it
//         # saves post-processing the (possibly considerable) expense of
//         # figuring out what to do with it.  In the case of an empty
//         # interesting match, this is clearly the right thing to do,
//         # because no other kind of match is possible in the regions.
//         while besti > alo and bestj > blo and \
//               isbjunk(b[bestj-1]) and \
//               a[besti-1] == b[bestj-1]:
//             besti, bestj, bestsize = besti-1, bestj-1, bestsize+1
//         while besti+bestsize < ahi and bestj+bestsize < bhi and \
//               isbjunk(b[bestj+bestsize]) and \
//               a[besti+bestsize] == b[bestj+bestsize]:
//             bestsize = bestsize + 1

//         return Match(besti, bestj, bestsize)

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