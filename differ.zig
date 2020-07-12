const std = @import("std");
const seq_matcher = @import("./sequence_matcher.zig");
const SequenceMatcher = seq_matcher.SequenceMatcher;
const OpcodeTag = seq_matcher.OpcodeTag;
const Opcode = seq_matcher.Opcode;
const mem = std.mem;
const fmt = std.fmt;
const testing = std.testing;
const assert = std.debug.assert;
const ArrayList = std.ArrayList;

const SetSeq2ErrorSet = @TypeOf(seq_matcher.SequenceMatcher([]const u8).setSeq2).ReturnType.ErrorSet;

pub const LineDiff = struct {
    from_line:    []const u8,
    to_line: []const u8,
    change: bool,
};

const LinesAndStarters = struct {
    lines: []const []const u8,
    starters: []const u8,
};

// Differ takes [][]const u8 lines of input (each should end with a \n)
// and outputs lines of Diffs ([][]const u8).
pub const Differ = struct {
    allocator: *mem.Allocator,
    sm: SequenceMatcher([]const []const u8),
    diffs: ArrayList([]u8),
    diff_lines: ArrayList(LineDiff),

    pub fn init(allocator: *mem.Allocator) Differ {
        return Differ{
            .allocator = allocator,
            .sm = SequenceMatcher([]const []const u8).init(allocator),
            .diffs = ArrayList([]u8).init(allocator),
            .diff_lines = ArrayList(LineDiff).init(allocator),
        };
    }

    pub fn deinit(self: *Differ) void {
        for (self.diffs.items) |diff| {
            self.allocator.free(diff);
        }
        self.diffs.deinit();
        self.sm.deinit();
    }

    pub fn reset(self: *Differ) void {
        for (self.diffs.items) |diff| {
            self.allocator.free(diff);
        }
        self.diffs.items.len = 0;
    }

    pub fn compare(self: *Differ, a: []const []const u8, b: []const []const u8) ![][]u8 {
        try self.sm.setSeqs(a, b);
        const opcodes = try self.sm.getOpcodes();
        for (opcodes) |opcode| {
            switch (opcode.tag) {
                OpcodeTag.Replace => {
                    try self.fancyReplace(a, b, opcode.a_start, opcode.a_end, opcode.b_start, opcode.b_end);
                },
                OpcodeTag.Delete => {
                    try self.dump("-", a, opcode.a_start, opcode.a_end);
                },
                OpcodeTag.Insert => {
                    try self.dump("+", b, opcode.b_start, opcode.b_end);
                },
                OpcodeTag.Equal => {
                    try self.dump(" ", a, opcode.a_start, opcode.a_end);
                },
            }
        }
        return self.diffs.items;
    }

    fn dump(self: *Differ, tag: []const u8, seq: []const []const u8, start: usize, end: usize) !void {
        var i: usize = start;
        while (i < end) : (i += 1) {
            var res = try fmt.allocPrint(self.allocator, "{} {}", .{ tag, seq[i] });
            try self.diffs.append(res);
        }
    }

    fn fancyReplace(self: *Differ, a: []const []const u8, b: []const []const u8, a_start: usize, a_end: usize, b_start: usize, b_end: usize) SetSeq2ErrorSet!void {
        var best_ratio: f32 = 0.74;
        var cutoff: f32 = 0.75;
        var best_i: usize = undefined;
        var best_j: usize = undefined;

        var sm = SequenceMatcher([]const u8).init(self.allocator);
        defer sm.deinit();

        var eqi: ?usize = null;
        var eqj: ?usize = null;

        var j: usize = b_start;
        while (j < b_end) : (j += 1) {
            var bj = b[j];
            try sm.setSeq2(bj);

            var i: usize = a_start;
            while (i < a_end) : (i += 1) {
                var ai = a[i];
                if (mem.eql(u8, ai, bj)) {
                    if (eqi == null) {
                        eqi = i;
                        eqj = j;
                    }
                    continue;
                }
                sm.setSeq1(ai);
                if (sm.realQuickRatio() > best_ratio and try sm.quickRatio() > best_ratio and try sm.ratio() > best_ratio) {
                    best_ratio = try sm.ratio();
                    best_i = i;
                    best_j = j;
                }
            }
        }

        if (best_ratio < cutoff) {
            if (eqi) |e| {
                best_i = e;
                best_j = eqj.?;
                best_ratio = 1.0;
            } else {
                try self.plainReplace(a, b, a_start, a_end, b_start, b_end);
                return;
            }
        } else {
            eqi = null;
        }


        if (a_start < best_i) {
            if (b_start < best_j) {
                try self.fancyReplace(a, b, a_start, best_i, b_start, best_j);
            } else {
                try self.dump("-", a, a_start, best_i);
            }
        } else if (b_start < best_j) {
            try self.dump("+", b, b_start, best_j);
        }

        var a_el = a[best_i];
        var b_el = b[best_j];
        if (eqi) |_| {
            var res = try fmt.allocPrint(self.allocator, "  {}", .{a_el});
            try self.diffs.append(res);
        } else {
            var a_tags = ArrayList(u8).init(self.allocator);
            var b_tags = ArrayList(u8).init(self.allocator);
            defer {
                a_tags.deinit();
                b_tags.deinit();
            }
            try sm.setSeqs(a_el, b_el);
            var opcodes = try sm.getOpcodes();
            for (opcodes) |opcode| {
                var la: usize = opcode.a_end - opcode.a_start;
                var lb: usize = opcode.b_end - opcode.b_start;
                switch (opcode.tag) {
                    OpcodeTag.Replace => {
                        var idx: usize = 0;
                        while (idx < la) : (idx += 1) {
                            try a_tags.append('^');
                        }
                        idx = 0;
                        while (idx < lb) : (idx += 1) {
                            try b_tags.append('^');
                        }
                    },
                    OpcodeTag.Delete => {
                        var idx: usize = 0;
                        while (idx < la) : (idx += 1) {
                            try a_tags.append('-');
                        }
                    },
                    OpcodeTag.Insert => {
                        var idx: usize = 0;
                        while (idx < lb) : (idx += 1) {
                            try b_tags.append('+');
                        }
                    },
                    OpcodeTag.Equal => {
                        var idx: usize = 0;
                        while (idx < la) : (idx += 1) {
                            try a_tags.append(' ');
                        }
                        idx = 0;
                        while (idx < lb) : (idx += 1) {
                            try b_tags.append(' ');
                        }
                    },
                }
            }
            try self.qformat(a_el, b_el, a_tags.items, b_tags.items);
        }
        if (best_i + 1 < a_end) {
            if (best_j + 1 < b_end) {
                try self.fancyReplace(a, b, best_i + 1, a_end, best_j + 1, b_end);
            } else {
                try self.dump("-", a, best_i + 1, a_end);
            }
        } else if (best_j + 1 < b_end) {
            try self.dump("+", b, best_j + 1, b_end);
        }
    }

    fn plainReplace(self: *Differ, a: []const []const u8, b: []const []const u8, a_start: usize, a_end: usize, b_start: usize, b_end: usize) !void {
        assert(a_start < a_end and b_start < b_end);

        if (b_end - b_start < a_end - a_start) {
            try self.dump("+", b, b_start, b_end);
            try self.dump("-", a, a_start, a_end);
        } else {
            try self.dump("-", a, a_start, a_end);
            try self.dump("+", b, b_start, b_end);
        }
    }

    fn qformat(self: *Differ, a_line: []const u8, b_line: []const u8, a_tags: []const u8, b_tags: []const u8) !void {
        var clean_a_tags = try keepOriginalWhitespace(self.allocator, a_line, a_tags);
        var clean_b_tags = try keepOriginalWhitespace(self.allocator, b_line, b_tags);
        defer {
            clean_a_tags.deinit();
            clean_b_tags.deinit();
        }

        var al = try fmt.allocPrint(self.allocator, "- {}", .{a_line});
        try self.diffs.append(al);

        if (clean_a_tags.items.len > 0) {
            var at = try fmt.allocPrint(self.allocator, "? {}\n", .{clean_a_tags.items});
            try self.diffs.append(at);
        }

        var bl = try fmt.allocPrint(self.allocator, "+ {}", .{b_line});
        try self.diffs.append(bl);

        if (clean_b_tags.items.len > 0) {
            var bt = try fmt.allocPrint(self.allocator, "? {}\n", .{clean_b_tags.items});
            try self.diffs.append(bt);
        }
    }

    fn getNextLines(self: *Differ, idx: usize) !LinesAndStarters {
        assert(idx < self.diffs.items.len);
        var lines = try ArrayList([]u8).initCapacity(self.allocator, 4);

        var x = [_]u8{'X'};

        var diffs = self.diffs.items;
        var i = idx;
        while (i < idx + 4) : (i += 1) {
            if (i >= diffs.len) {
                try lines.append(&x);
            } else {
                try lines.append(diffs[i]);
            }
        }

        var starters = try ArrayList(u8).initCapacity(self.allocator, 4);
        for (lines.items) |line| {
            try starters.append(line[0]);
        }

        return LinesAndStarters{
            .lines = lines.toOwnedSlice(),
            .starters = starters.toOwnedSlice(),
        };
    }

    fn getDiffLines(self: *Differ) !void {
        assert(self.diffs.items.len > 0);

        var blanks_pending: i32 = 0;
        var blanks_ready: i32 = 0;

        var diff_idx: usize = 0;
        while (diff_idx < self.diffs.items.len) : (diff_idx += 4) {
            var lines_and_starters = try self.getNextLines(diff_idx);
            defer {
                self.allocator.free(lines_and_starters.starters);
                self.allocator.free(lines_and_starters.lines);
            }
            var s = lines_and_starters.starters;

            if (mem.startsWith(u8, s, "X")) {
                blanks_ready = blanks_pending;
            } else if (mem.startsWith(u8, s, "-?+?")) {
                // make line
                continue;
            } else if (mem.startsWith(u8, s, "--++")) {
                blanks_pending -= 1;
                // make line
                continue;
            } else if (mem.startsWith(u8, s, "--?+") or mem.startsWith(u8, s, "--+") or mem.startsWith(u8, s, "- ")) {
                // make line
                blanks_ready = blanks_pending - 1;
                blanks_pending = 0;
            } else if (mem.startsWith(u8, s, "-+?")) {
                // make line
                continue;
            } else if (mem.startsWith(u8, s, "-?+")) {
                // make line
                continue;
            } else if (mem.startsWith(u8, s, "-")) {
                blanks_pending -= 1;
                // make line
                continue;
            } else if (mem.startsWith(u8, s, "+--")) {
                blanks_pending += 1;
                // make line
                continue;
            } else if (mem.startsWith(u8, s, "+ ") or mem.startsWith(u8, s, "+-")) {
                // make line
                blanks_ready = blanks_pending + 1;
                blanks_pending = 0;
            } else if (mem.startsWith(u8, s, "+")) {
                blanks_pending += 1;
                // make line
                continue;
            } else if (mem.startsWith(u8, s, " ")) {
                // make line
                continue;
            }

            while (blanks_ready < 0) : (blanks_ready += 1) {
                // yield something
            }

            while (blanks_pending > 0) : (blanks_pending -= 1) {
                // yield something
            }

            if (mem.startsWith(u8, s, "X")) {
                return; // return what?
            } else {
                // yield something
            }

        }
    }

    pub fn mdiff(self: *Differ, from_lines: []const []const u8, to_lines: []const []const u8) !void {
        var diffs = try self.compare(from_lines, to_lines);
        try self.getDiffLines();
    }
};

fn keepOriginalWhitespace(allocator: *mem.Allocator, str: []const u8, tag_str: []const u8) !ArrayList(u8) {
    var out = ArrayList(u8).init(allocator);
    for (str) |c, idx| {
        if (idx >= tag_str.len) break;
        var tag_c = tag_str[idx];
        if (tag_c == ' ' and std.ascii.isSpace(c)) {
            try out.append(c);
        } else {
            try out.append(tag_c);
        }
    }
    // remove right-hand side whitespace
    var idx: usize = out.items.len;
    while (idx > 0 and std.ascii.isSpace(out.items[idx - 1])) {
        idx -= 1;
    }
    out.items.len = idx;
    return out;
}

test "mdiff" {
    var allocator = testing.allocator;

    var differ = Differ.init(allocator);
    defer differ.deinit();

    const a = &[_][]const u8{ "one\n", "two\n", "three\n" };
    const b = &[_][]const u8{ "ore\n", "tree\n", "emu\n" };
    try differ.mdiff(a, b);
}


test "compare" {
    const TestCase = struct {
        a: []const []const u8,
        b: []const []const u8,
        expected: []const []const u8,
    };

    const testCases = [_]TestCase{
        TestCase{
            .a = &[_][]const u8{ "one\n", "two\n", "three\n" },
            .b = &[_][]const u8{ "ore\n", "tree\n", "emu\n" },
            .expected = &[_][]const u8{ "- one\n", "?  ^\n", "+ ore\n", "?  ^\n", "- two\n", "- three\n", "?  -\n", "+ tree\n", "+ emu\n" },
        },
        TestCase{
            .a = &[_][]const u8{ "\tI am a buggy" },
            .b = &[_][]const u8{ "\t\tI am a bug" },
            .expected = &[_][]const u8{ "- \tI am a buggy", "? \t          --\n", "+ \t\tI am a bug", "? +\n" },
        },
        TestCase{
            .a = &[_][]const u8{ "\t \t \t^" },
            .b = &[_][]const u8{ "\t \t \t^\n" },
            .expected = &[_][]const u8{ "- \t \t \t^", "+ \t \t \t^\n", "? \t \t \t +\n" },
        },
    };

    const allocator = testing.allocator;

    var differ = Differ.init(allocator);
    defer differ.deinit();

    for (testCases) |tc| {
        var diffs = try differ.compare(tc.a[0..], tc.b[0..]);
        defer differ.reset();
        assert(diffs.len == tc.expected.len);
        for (tc.expected) |exp_line, idx| {
            var actual_line = diffs[idx];
            assert(mem.eql(u8, exp_line, actual_line));
        }
    }
}
