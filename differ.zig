const std = @import("std");
const seq_matcher = @import("./sequence_matcher.zig");
const SequenceMatcher = seq_matcher.SequenceMatcher;
const OpcodeTag = seq_matcher.OpcodeTag;
const Opcode = seq_matcher.Opcode;
const mem = std.mem;
const fmt = std.fmt;
const testing = std.testing;
const ArrayList = std.ArrayList;

pub const Differ = struct {
    allocator: *mem.Allocator,
    sm: SequenceMatcher,
    opcodes: ?[]Opcode,
    opcodes_idx: usize,
    frame: ?anyframe,
    next: ?[]u8,

    pub fn init(allocator: *mem.Allocator) !Differ {
        return Differ{
            .allocator = allocator,
            .sm = try SequenceMatcher.init(allocator, "", ""),
            .opcodes = null,
            .opcodes_idx = 0,
            .frame = null,
            .next = null,
        };
    }

    pub fn deinit(self: *Differ) void {
        self.sm.deinit();
    }

    pub const CompareGenerator = struct {
        allocator: *mem.Allocator,
        a: []const u8,
        b: []const u8,
        opcodes: []Opcode,
        curr_opcode: usize = 0,
        dump_gen: ?DumpGenerator = null,
        repl_gen: ?DumpGenerator = null,

        pub fn next(gen: *CompareGenerator) !?[]u8 {
            var res = try gen.takeFromGenerator();
            if (res) |_| {
                return res;
            } else {
                // our internal generators are exhausted;
                // reset everything and move to next opcode
                gen.dump_gen = null;
                gen.repl_gen = null;
                gen.initGenerator();
                var _res = try gen.takeFromGenerator();
                gen.curr_opcode += 1;
                return _res;
            }
        }

        fn takeFromGenerator(gen: *CompareGenerator) !?[]u8 {
            var res: ?[]u8 = null;
            if (gen.dump_gen) |*dump_gen| {
                res = try dump_gen.next();
            } else if (gen.repl_gen) |*repl_gen| {
                res = try repl_gen.next();
            }
            return res;
        }

        fn initGenerator(gen: *CompareGenerator) void {
            if (gen.curr_opcode >= gen.opcodes.len) return;
            var opcode = gen.opcodes[gen.curr_opcode];

            var next_dump_gen: ?DumpGenerator = null;
            var next_repl_gen: ?DumpGenerator = null;

            switch (opcode.tag) {
                OpcodeTag.Delete => {
                    gen.dump_gen = DumpGenerator{
                        .start = opcode.a_start,
                        .end = opcode.a_end,
                        .tag = "-",
                        .seq = gen.a,
                        .allocator = gen.allocator,
                    };
                },
                OpcodeTag.Insert => {
                    gen.dump_gen = DumpGenerator{
                        .start = opcode.b_start,
                        .end = opcode.b_end,
                        .tag = "+",
                        .seq = gen.b,
                        .allocator = gen.allocator,
                    };
                },
                OpcodeTag.Equal => {
                    gen.dump_gen = DumpGenerator{
                        .start = opcode.a_start,
                        .end = opcode.a_end,
                        .tag = " ",
                        .seq = gen.a,
                        .allocator = gen.allocator,
                    };
                },
                else => {},
            }
        }
    };

    pub fn compare(self: *Differ, a: []const u8, b: []const u8) !CompareGenerator {
        try self.sm.setSeqs(a, b);
        var opcodes = try self.sm.getOpcodes();
        return CompareGenerator{
            .a = a,
            .b = b,
            .opcodes = opcodes,
            .allocator = self.allocator,
        };
    }

    pub const DumpGenerator = struct {
        start: usize,
        end: usize,
        tag: []const u8,
        seq: []const u8,
        allocator: *mem.Allocator,

        pub fn next(gen: *DumpGenerator) !?[]u8 {
            if (gen.start >= gen.end) return null;
            var res = try fmt.allocPrint(gen.allocator, "{} {}", .{gen.tag, gen.seq[gen.start..gen.start+1]});
            gen.start += 1;

            return res;
        }
    };
};

fn collectLines(allocator: *mem.Allocator, it: var) !ArrayList([]u8) {
    var output = ArrayList([]u8).init(allocator);
    while (try it.next()) |item| {
        try output.append(item);
    }
    return output;
}

test "compare" {
    var differ = try Differ.init(testing.allocator);
    defer differ.deinit();

    var it = try differ.compare("asdf\nfdsa\n", "asxdf\na\n");
    while (try it.next()) |diff_line| {
        std.debug.warn("\n{}", .{diff_line});
        testing.allocator.free(diff_line);
    }
    std.debug.warn("\n", .{});
}
