const std = @import("std");
const mem = std.mem;
const assert = std.debug.assert;
const sort = std.sort;
const testing = std.testing;
const ArrayList = std.ArrayList;
const TailQueue = std.TailQueue;
const AutoHashMap = std.AutoHashMap;
const StringHashMap = std.StringHashMap;

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

pub const OpcodeTag = enum {
    Replace,
    Delete,
    Insert,
    Equal,
};

pub const Opcode = struct {
    tag: OpcodeTag,
    a_start: usize,
    a_end: usize,
    b_start: usize,
    b_end: usize,
};

fn calculateRatio(matches: usize, length: usize) f32 {
    var m = @intToFloat(f32, matches);
    var l = @intToFloat(f32, length);
    if (length > 0) {
        return 2.0 * m / l;
    }
    return 1.0;
}

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

pub fn SequenceMatcher(comptime T: type) type {
    assert(T == []const u8 or T == []const []const u8);

    const SeqType = T;
    const ElType = if (T == []const u8) u8 else []const u8;
    const SeqToIdxMapType = if (T == []const u8) AutoHashMap(ElType, ArrayList(usize)) else StringHashMap(ArrayList(usize));
    const SeqToCountMapType = if (T == []const u8) AutoHashMap(ElType, i32) else StringHashMap(i32);

    return struct {
        const Self = @This();
        seq1: SeqType,
        seq2: SeqType,
        allocator: *mem.Allocator,

        b2j: SeqToIdxMapType,
        matching_blocks: ?ArrayList(Match),
        opcodes: ?ArrayList(Opcode),
        full_b_count: ?SeqToCountMapType,

        pub fn init(allocator: *mem.Allocator) Self {
            return Self{
                .allocator = allocator,
                .seq1 = &[0]ElType{},
                .seq2 = &[0]ElType{},
                .b2j = SeqToIdxMapType.init(allocator),
                .matching_blocks = null,
                .opcodes = null,
                .full_b_count = null,
            };
        }

        pub fn deinit(self: *Self) void {
            var it = self.b2j.iterator();
            while (it.next()) |indices| {
                indices.value.deinit();
            }
            self.b2j.deinit();
            if (self.matching_blocks) |mb| {
                mb.deinit();
            }
            if (self.opcodes) |opcodes| {
                opcodes.deinit();
            }
            if (self.full_b_count) |fbc| {
                fbc.deinit();
            }
        }

        pub fn setSeqs(self: *Self, a: SeqType, b: SeqType) !void {
            self.setSeq1(a);
            try self.setSeq2(b);
        }

        pub fn setSeq1(self: *Self, a: SeqType) void {
            if (self.equalSeqs(a, self.seq1)) {
                return;
            }

            self.seq1 = a;
            if (self.matching_blocks) |mb| {
                mb.deinit();
            }
            self.matching_blocks = null;
            if (self.opcodes) |opcodes| {
                opcodes.deinit();
            }
            self.opcodes = null;
        }

        pub fn setSeq2(self: *Self, b: SeqType) !void {
            if (self.equalSeqs(b, self.seq2)) {
                return;
            }
            self.seq2 = b;
            if (self.matching_blocks) |mb| {
                mb.deinit();
            }
            self.matching_blocks = null;
            if (self.opcodes) |opcodes| {
                opcodes.deinit();
            }
            self.opcodes = null;
            if (self.full_b_count) |fbc| {
                fbc.deinit();
            }
            self.full_b_count = null;

            try self.chainB();
        }

        fn clearB2j(self: *Self) void {
            var it = self.b2j.iterator();
            while (it.next()) |indices| {
                indices.value.deinit();
            }
            self.b2j.clear();
        }

        fn chainB(self: *Self) !void {
            var b = self.seq2;
            self.clearB2j();

            for (b) |el, idx| {
                if (self.b2j.getValue(el)) |*idxs| {
                    try idxs.append(idx);
                    _ = try self.b2j.put(el, idxs.*);
                } else {
                    var idxs = ArrayList(usize).init(self.allocator);
                    try idxs.append(idx);
                    _ = try self.b2j.put(el, idxs);
                }
            }
            if (b.len >= 200) {
                var popular_els = ArrayList(ElType).init(self.allocator);
                defer popular_els.deinit();
                var n_test = @intToFloat(f32, b.len) / 100 + 1;
                var it = self.b2j.iterator();
                while (it.next()) |kv| {
                    var el_idx_size = @intToFloat(f32, kv.value.items.len);
                    if (el_idx_size > n_test) {
                        kv.value.deinit();
                        try popular_els.append(kv.key);
                    }
                }
                for (popular_els.items) |popular_el| {
                    _ = self.b2j.remove(popular_el);
                }
            }
        }

        fn findLongestMatch(self: *Self, a_lo: usize, a_hi: usize, b_lo: usize, b_hi: usize) !Match {
            var a = self.seq1;
            var b = self.seq2;
            assert(a_hi <= a.len and b_hi <= b.len);

            var b2j = self.b2j;
            var best_i = a_lo;
            var best_j = b_lo;
            var best_size: usize = 0;
            var j2_len = AutoHashMap(usize, usize).init(self.allocator);
            var j2_len_next = AutoHashMap(usize, usize).init(self.allocator);
            defer {
                j2_len.deinit();
                j2_len_next.deinit();
            }
            var nothing = [0]usize{};

            var i = a_lo;
            while (i < a_hi) : (i += 1) {
                var indices: []usize = undefined;
                if (b2j.getValue(a[i])) |value| {
                    indices = value.items;
                } else {
                    indices = &nothing;
                }
                for (indices) |j| {
                    if (j < b_lo) continue;
                    if (j >= b_hi) break;
                    var k = (if (j == 0) 0 else (j2_len.getValue(j - 1) orelse 0)) + 1;
                    _ = try j2_len_next.put(j, k);
                    if (k > best_size) {
                        best_i = i + 1 - k;
                        best_j = j + 1 - k;
                        best_size = k;
                    }
                }
                var tmp = j2_len;
                j2_len = j2_len_next;
                j2_len_next = tmp;
                j2_len_next.clear();
            }

            while (best_i > a_lo and best_j > b_lo and self.equalEls(a[best_i - 1], b[best_j - 1])) {
                best_i -= 1;
                best_j -= 1;
                best_size += 1;
            }

            while (best_i + best_size < a_hi and best_j + best_size < b_hi and self.equalEls(a[best_i + best_size], b[best_j + best_size])) {
                best_size += 1;
            }

            return Match{
                .a_idx = best_i,
                .b_idx = best_j,
                .size = best_size,
            };
        }

        pub fn getMatchingBlocks(self: *Self) ![]Match {
            if (self.matching_blocks) |mb| {
                return mb.items;
            }

            var len_a = self.seq1.len;
            var len_b = self.seq2.len;

            var matching_blocks = ArrayList(Match).init(self.allocator);
            defer matching_blocks.deinit();
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

            var a_idx: usize = 0;
            var b_idx: usize = 0;
            var size: usize = 0;
            var non_adjacent = ArrayList(Match).init(self.allocator);
            errdefer non_adjacent.deinit();
            for (matching_blocks.items) |match| {
                if (a_idx + size == match.a_idx and b_idx + size == match.b_idx) {
                    size += match.size;
                } else {
                    if (size > 0) {
                        try non_adjacent.append(Match{
                            .a_idx = a_idx,
                            .b_idx = b_idx,
                            .size = size,
                        });
                    }
                    a_idx = match.a_idx;
                    b_idx = match.b_idx;
                    size = match.size;
                }
            }
            if (size > 0) {
                try non_adjacent.append(Match{
                    .a_idx = a_idx,
                    .b_idx = b_idx,
                    .size = size,
                });
            }

            // sentinal match
            try non_adjacent.append(Match{
                .a_idx = len_a,
                .b_idx = len_b,
                .size = 0,
            });

            self.matching_blocks = non_adjacent;
            return non_adjacent.items;
        }

        pub fn getOpcodes(self: *Self) ![]Opcode {
            if (self.opcodes) |opcodes| {
                return opcodes.items;
            }
            var opcodes = ArrayList(Opcode).init(self.allocator);

            var i: usize = 0;
            var j: usize = 0;

            var mbs = try self.getMatchingBlocks();
            for (mbs) |match| {
                var tag: ?OpcodeTag = null;
                if (i < match.a_idx and j < match.b_idx) {
                    tag = OpcodeTag.Replace;
                } else if (i < match.a_idx) {
                    tag = OpcodeTag.Delete;
                } else if (j < match.b_idx) {
                    tag = OpcodeTag.Insert;
                }
                if (tag) |t| {
                    try opcodes.append(Opcode{
                        .tag = t,
                        .a_start = i,
                        .a_end = match.a_idx,
                        .b_start = j,
                        .b_end = match.b_idx,
                    });
                }
                i = match.a_idx + match.size;
                j = match.b_idx + match.size;
                if (match.size > 0) {
                    try opcodes.append(Opcode{
                        .tag = OpcodeTag.Equal,
                        .a_start = match.a_idx,
                        .a_end = i,
                        .b_start = match.b_idx,
                        .b_end = j,
                    });
                }
            }

            self.opcodes = opcodes;
            return opcodes.items;
        }

        pub fn ratio(self: *Self) !f32 {
            var sum_of_matches: usize = 0;
            for (try self.getMatchingBlocks()) |match| {
                sum_of_matches += match.size;
            }
            return calculateRatio(sum_of_matches, self.seq1.len + self.seq2.len);
        }

        pub fn quickRatio(self: *Self) !f32 {
            var full_b_count: SeqToCountMapType = undefined;
            if (self.full_b_count) |fbc| {
                full_b_count = fbc;
            } else {
                full_b_count = SeqToCountMapType.init(self.allocator);
                for (self.seq2) |el| {
                    if (full_b_count.getValue(el)) |count| {
                        _ = try full_b_count.put(el, count + 1);
                    } else {
                        _ = try full_b_count.put(el, 1);
                    }
                }
                self.full_b_count = full_b_count;
            }

            var avail = SeqToCountMapType.init(self.allocator);
            defer avail.deinit();
            var matches: usize = 0;
            for (self.seq1) |el| {
                var n: i32 = 0;
                if (avail.getValue(el)) |a| {
                    n = a;
                } else {
                    n = full_b_count.getValue(el) orelse 0;
                }
                _ = try avail.put(el, n - 1);
                if (n > 0) matches += 1;
            }
            return calculateRatio(matches, self.seq1.len + self.seq2.len);
        }

        pub fn realQuickRatio(self: *Self) f32 {
            var len_a = self.seq1.len;
            var len_b = self.seq2.len;
            return calculateRatio(if (len_a < len_b) len_a else len_b, len_a + len_b);
        }

        fn equalSeqs(self: *Self, a: SeqType, b: SeqType) bool {
            if (SeqType == []const u8) {
                return mem.eql(u8, a, b);
            }
            if (a.len != b.len) return false;
            for (a) |el, idx| {
                if (!mem.eql(u8, el, b[idx])) return false;
            }
            return true;
        }

        fn equalEls(self: *Self, a: ElType, b: ElType) bool {
            if (ElType == u8) {
                return a == b;
            } else {
                return mem.eql(u8, a, b);
            }
        }
    };
}

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

    const allocator = testing.allocator;
    var sm = SequenceMatcher([]const u8).init(allocator);
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
    var allocator = testing.allocator;
    var sm = SequenceMatcher([]const u8).init(allocator);
    try sm.setSeqs("abxcd", "abcd");
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

test "get opcodes" {
    const TestCase = struct {
        a: []const u8,
        b: []const u8,
        expected: []const Opcode,
    };

    const testCases = [_]TestCase{
        TestCase{
            .a = "b" ** 100,
            .b = "a" ++ "b" ** 100,
            .expected = &[_]Opcode{
                Opcode{
                    .tag = OpcodeTag.Insert,
                    .a_start = 0,
                    .a_end = 0,
                    .b_start = 0,
                    .b_end = 1,
                },
                Opcode{
                    .tag = OpcodeTag.Equal,
                    .a_start = 0,
                    .a_end = 100,
                    .b_start = 1,
                    .b_end = 101,
                },
            },
        },
        TestCase{
            .a = "b" ** 100,
            .b = "b" ** 50 ++ "a" ++ "b" ** 50,
            .expected = &[_]Opcode{
                Opcode{
                    .tag = OpcodeTag.Equal,
                    .a_start = 0,
                    .a_end = 50,
                    .b_start = 0,
                    .b_end = 50,
                },
                Opcode{
                    .tag = OpcodeTag.Insert,
                    .a_start = 50,
                    .a_end = 50,
                    .b_start = 50,
                    .b_end = 51,
                },
                Opcode{
                    .tag = OpcodeTag.Equal,
                    .a_start = 50,
                    .a_end = 100,
                    .b_start = 51,
                    .b_end = 101,
                },
            },
        },
        TestCase{
            .a = "a" ** 40 ++ "c" ++ "b" ** 40,
            .b = "a" ** 40 ++ "b" ** 40,
            .expected = &[_]Opcode{
                Opcode{
                    .tag = OpcodeTag.Equal,
                    .a_start = 0,
                    .a_end = 40,
                    .b_start = 0,
                    .b_end = 40,
                },
                Opcode{
                    .tag = OpcodeTag.Delete,
                    .a_start = 40,
                    .a_end = 41,
                    .b_start = 40,
                    .b_end = 40,
                },
                Opcode{
                    .tag = OpcodeTag.Equal,
                    .a_start = 41,
                    .a_end = 81,
                    .b_start = 40,
                    .b_end = 80,
                },
            },
        },
    };

    var sm = SequenceMatcher([]const u8).init(testing.allocator);
    defer sm.deinit();
    for (testCases) |testCase| {
        try sm.setSeqs(testCase.a, testCase.b);
        var opcodes = try sm.getOpcodes();
        assert(opcodes.len == testCase.expected.len);
        for (opcodes) |opcode, idx| {
            var expected = testCase.expected[idx];
            assert(opcode.tag == expected.tag);
            assert(opcode.a_start == expected.a_start);
            assert(opcode.a_end == expected.a_end);
            assert(opcode.b_start == expected.b_start);
            assert(opcode.b_end == expected.b_end);
        }
    }
}

test "ratio" {
    const TestCase = struct {
        a: []const u8,
        b: []const u8,
        expected: f32,
    };

    const testCases = [_]TestCase{
        TestCase{
            .a = "b" ** 100,
            .b = "a" ++ "b" ** 100,
            .expected = 0.995,
        },
        TestCase{
            .a = "b" ** 100,
            .b = "b" ** 50 ++ "a" ++ "b" ** 50,
            .expected = 0.995,
        },
        TestCase{
            .a = "a" ** 40 ++ "c" ++ "b" ** 40,
            .b = "a" ** 40 ++ "b" ** 40,
            .expected = 0.994,
        },
        TestCase{
            .a = "b" ** 200,
            .b = "a" ++ "b" ** 200,
            .expected = 0.0,
        },
        TestCase{
            .a = "",
            .b = "",
            .expected = 1.0,
        },
    };

    var sm = SequenceMatcher([]const u8).init(testing.allocator);
    defer sm.deinit();
    for (testCases) |testCase| {
        try sm.setSeqs(testCase.a, testCase.b);
        var actual = try sm.ratio();
        var quick = try sm.quickRatio();
        var real_quick = sm.realQuickRatio();
        assert(std.math.approxEq(f32, actual, testCase.expected, 0.001));
    }
}

test "quick ratio" {
    var sm = SequenceMatcher([]const u8).init(testing.allocator);
    defer sm.deinit();

    try sm.setSeqs("one\n", "ore\n");
    var quick = try sm.quickRatio();

    assert(quick == 0.75);
}

test "sequence matching of strings" {
    var sm = SequenceMatcher([]const []const u8).init(testing.allocator);
    defer sm.deinit();

    const a = [_][]const u8{
        "one\n",
        "two\n",
        "three\n",
    };
    const b = [_][]const u8{
        "ore\n",
        "tree\n",
        "emu\n",
    };

    try sm.setSeqs(a[0..], b[0..]);
    const opcodes = try sm.getOpcodes();
    assert(opcodes.len == 1);
    assert(opcodes[0].a_start == 0);
    assert(opcodes[0].a_end == 3);
    assert(opcodes[0].b_start == 0);
    assert(opcodes[0].b_end == 3);
}
