//! Builds libskar.{a,lib}: a static archive that exposes the skar Zig
//! package's `solve` via a C ABI for the Cython extension to link
//! against. The upstream skar source is resolved from the dependency
//! pinned in build.zig.zon — a local `.path` during development, a
//! URL+hash for releases.

const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const skar_mod = b.dependency("skar", .{
        .target = target,
        .optimize = optimize,
    }).module("skar");

    const cabi_mod = b.createModule(.{
        .root_source_file = b.path("c_api.zig"),
        .target = target,
        .optimize = optimize,
        // skar.solve allocates; the shim hands it std.heap.c_allocator.
        .link_libc = true,
        // The static archive ends up linked into a Python extension
        // (.so / .pyd), itself a shared library — its objects must be
        // position-independent.
        .pic = true,
        .imports = &.{
            .{ .name = "skar", .module = skar_mod },
        },
    });

    // Static lib pulled into the Cython extension at link time. Avoids
    // the Windows MSVC CRT mismatch and the macOS dylib __dso_handle
    // regression that shipping a Zig *dynamic* library triggers — the
    // same rationale documented in the sibling sparea_py bindings.
    const lib = b.addLibrary(.{
        .name = "skar",
        .linkage = .static,
        .root_module = cabi_mod,
    });
    b.installArtifact(lib);
}
