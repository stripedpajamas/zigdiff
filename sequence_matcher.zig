const std = @import("std");
const mem = std.mem;
const assert = std.debug.assert;
const sort = std.sort;
const ArrayList = std.ArrayList;
const TailQueue = std.TailQueue;
const AutoHashMap = std.AutoHashMap;

pub const Match = struct {
    a_idx: usize,
    b_idx: usize,
    size: usize,
};
const MatchParams = struct {
    a_lo: usize,
    a_hi: usize,
    b_lo: usize,
    b_hi: usize,
};

fn matchLessThan(match_a: Match, match_b: Match) bool {
    // is match_a < match_b ?
    if (match_a.a_idx < match_b.a_idx) return true;
    if (match_a.a_idx > match_b.a_idx) return false;
    if (match_a.a_idx == match_b.a_idx) {
        if (match_a.b_idx < match_b.b_idx) return true;
        if (match_a.b_idx > match_b.b_idx) return false;
    }
    return true;
}

pub const SequenceMatcher = struct {
    seq1: []const u8,
    seq2: []const u8,
    allocator: *mem.Allocator,

    b2j: AutoHashMap(u8, ArrayList(usize)),
    matching_blocks: ?ArrayList(Match),
    opcodes: []i32,
    full_b_count: []i32,

    pub fn init(allocator: *mem.Allocator, a: []const u8, b: []const u8) !SequenceMatcher {
        var b2j = AutoHashMap(u8, ArrayList(usize)).init(allocator);
        var sm = SequenceMatcher{
            .allocator = allocator,
            .seq1 = undefined,
            .seq2 = undefined,
            .b2j = b2j,
            .matching_blocks = null,
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
        if (self.matching_blocks) |mb| {
            mb.deinit();
        }
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
        if (self.matching_blocks) |mb| {
            mb.deinit();
        }
        self.matching_blocks = null;
        self.opcodes = undefined;
    }

    pub fn setSeq2(self: *SequenceMatcher, b: []const u8) !void {
        if (mem.eql(u8, b, self.seq2)) {
            return;
        }
        self.seq2 = b;
        if (self.matching_blocks) |mb| {
            mb.deinit();
        }
        self.matching_blocks = null;
        self.opcodes = undefined;
        self.full_b_count = undefined;

        try self.chainB();
    }

    fn clearB2j(self: *SequenceMatcher) void {
        var it = self.b2j.iterator();
        while (it.next()) |indices| {
            indices.value.deinit();
        }
        self.b2j.clear();
    }

    fn chainB(self: *SequenceMatcher) !void {
        var b = self.seq2;
        self.clearB2j();

        for (b) |byte, idx| {
            if (self.b2j.getValue(byte)) |*byte_idxs| {
                try byte_idxs.append(idx);
                _ = try self.b2j.put(byte, byte_idxs.*);
            } else {
                var byte_idxs = ArrayList(usize).init(self.allocator);
                try byte_idxs.append(idx);
                _ = try self.b2j.put(byte, byte_idxs);
            }
        }
        if (b.len >= 200) {
            var popular_bytes = ArrayList(u8).init(self.allocator);
            defer popular_bytes.deinit();
            var n_test = @intToFloat(f32, b.len) / 100 + 1;
            var it = self.b2j.iterator();
            while (it.next()) |kv| {
                var byte_idx_size = @intToFloat(f32, kv.value.items.len);
                if (byte_idx_size > n_test) {
                    kv.value.deinit();
                    try popular_bytes.append(kv.key);
                }
            }
            for (popular_bytes.items) |popular_byte| {
                _ = self.b2j.remove(popular_byte);
            }
        }
    }

    fn findLongestMatch(self: *SequenceMatcher, a_lo: usize, a_hi: usize, b_lo: usize, b_hi: usize) !Match {
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
            j2_len.clear();
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

        return Match{
            .a_idx = best_i,
            .b_idx = best_j,
            .size = best_size,
        };
    }

    pub fn getMatchingBlocks(self: *SequenceMatcher) ![]Match {
        if (self.matching_blocks) |mb| {
            return mb.items;
        }
        
        var len_a = self.seq1.len;
        var len_b = self.seq2.len;

        var matching_blocks = ArrayList(Match).init(self.allocator);
        var queue = TailQueue(MatchParams).init();
        var initial_params = try queue.createNode(MatchParams{
            .a_lo = 0,
            .a_hi = len_a,
            .b_lo = 0,
            .b_hi = len_b,
        }, self.allocator);
        queue.append(initial_params);

        var prev_params = initial_params;
        while (queue.pop()) |queue_node| {
            var params = queue_node.data;
            var match = try self.findLongestMatch(params.a_lo, params.a_hi, params.b_lo, params.b_hi);
            if (match.size > 0) {
                try matching_blocks.append(match);
                if (params.a_lo < match.a_idx and params.b_lo < match.b_idx) {
                    var next_params = try queue.createNode(MatchParams{
                        .a_lo = params.a_lo,
                        .a_hi = match.a_idx,
                        .b_lo = params.b_lo,
                        .b_hi = match.b_idx,
                    }, self.allocator);
                    queue.append(next_params);
                }
                if (match.a_idx + match.size < params.a_hi and match.b_idx + match.size < params.b_hi) {
                    var next_params = try queue.createNode(MatchParams{
                        .a_lo = match.a_idx + match.size,
                        .a_hi = params.a_hi,
                        .b_lo = match.b_idx + match.size,
                        .b_hi = params.b_hi,
                    }, self.allocator);
                    queue.append(next_params);
                }
            }
            queue.destroyNode(prev_params, self.allocator);
            prev_params = queue_node;
        }

        sort.sort(Match, matching_blocks.items, matchLessThan);

        try matching_blocks.append(Match{
            .a_idx = len_a,
            .b_idx = len_b,
            .size = 0,
        });

        self.matching_blocks = matching_blocks;
        return matching_blocks.items;
    }

};

test "find longest match" {
    const TestCase = struct {
        a: []const u8,
        b: []const u8,
        a_lo: usize,
        a_hi: usize,
        b_lo: usize,
        b_hi: usize,
        expected: Match,
    };
    const testCases = [_]TestCase{
        TestCase{
            .a = "foo bar",
            .b = "foo baz bar",
            .a_lo = 0,
            .a_hi = 7,
            .b_lo = 0,
            .b_hi = 11,
            .expected = Match{
                .a_idx = 0,
                .b_idx = 0,
                .size = 6,
            },
        },
        TestCase{
            .a = "foo bar",
            .b = "foo baz bar",
            .a_lo = 2,
            .a_hi = 7,
            .b_lo = 4,
            .b_hi = 11,
            .expected = Match{
                .a_idx = 3,
                .b_idx = 7,
                .size = 4,
            },
        },
        TestCase{
            .a = "foo bar",
            .b = "foo baz bar",
            .a_lo = 0,
            .a_hi = 7,
            .b_lo = 1,
            .b_hi = 5,
            .expected = Match{
                .a_idx = 1,
                .b_idx = 1,
                .size = 4,
            },
        },
        TestCase{
            .a = "dabcd",
            .b = [_]u8{'d'} ** 100 ++ "abc" ++ [_]u8{'d'} ** 100,
            .a_lo = 0,
            .a_hi = 5,
            .b_lo = 0,
            .b_hi = 203,
            .expected = Match{
                .a_idx = 0,
                .b_idx = 99,
                .size = 5,
            },
        },
        TestCase{
            .a = "abcde",
            .b = "fghijk",
            .a_lo = 0,
            .a_hi = 5,
            .b_lo = 0,
            .b_hi = 6,
            .expected = Match{
                .a_idx = 0,
                .b_idx = 0,
                .size = 0,
            },
        },
        TestCase{
            .a = "abxcd",
            .b = "abcd",
            .a_lo = 0,
            .a_hi = 5,
            .b_lo = 0,
            .b_hi = 4,
            .expected = Match{
                .a_idx = 0,
                .b_idx = 0,
                .size = 2,
            },
        },
    };

    const allocator = std.testing.allocator;
    var sm = try SequenceMatcher.init(allocator, "", "");
    defer sm.deinit();

    for (testCases) |testCase| {
        try sm.setSeqs(testCase.a, testCase.b);
        var lm = try sm.findLongestMatch(testCase.a_lo, testCase.a_hi, testCase.b_lo, testCase.b_hi);
        assert(lm.a_idx == testCase.expected.a_idx);
        assert(lm.b_idx == testCase.expected.b_idx);
        assert(lm.size == testCase.expected.size);
    }
}

test "get matching blocks" {
    var allocator = std.testing.allocator;
    var sm = try SequenceMatcher.init(allocator, "abxcd", "abcd");
    defer sm.deinit();
    
    var mbs = try sm.getMatchingBlocks();
    var expected = [_]Match{
        Match{ .a_idx = 0, .b_idx = 0, .size = 2 },
        Match{ .a_idx = 3, .b_idx = 2, .size = 2 },
        Match{ .a_idx = 5, .b_idx = 4, .size = 0 },
    };

    assert(expected.len == mbs.len);
    for (mbs) |match, idx| {
        assert(match.a_idx == expected[idx].a_idx);
        assert(match.b_idx == expected[idx].b_idx);
        assert(match.size == expected[idx].size);
    }
}