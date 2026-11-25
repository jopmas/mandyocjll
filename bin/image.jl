# Compiled image of mandyocjll (built by `make`), shared by bin/mandyocjll and build/build_image.jl.
const ROOT = dirname(dirname(realpath(@__FILE__)))
const IMAGE = joinpath(ROOT, "bin", "mandyocjll." * (Sys.isapple() ? "dylib" : "so"))
const STAMP = joinpath(ROOT, "bin", "mandyocjll.stamp")   # 1st line: Julia that built the image; 2nd: checksum

"Changes when the code, the packages or the PETSc/MPI preferences change (independent of the Julia version)."
function image_checksum()
    files = [filter(f -> endswith(f, ".jl"), readdir(joinpath(ROOT, "src"); join = true));
        joinpath.(ROOT, ["Project.toml", "Manifest.toml", "LocalPreferences.toml"])]
    crc = UInt32(0)
    for f in files
        isfile(f) && (crc = Base._crc32c([codeunits(basename(f)); read(f)], crc))
    end
    return string(crc; base = 16)
end

"Julia executable that can use the image, or nothing."
function image_julia()
    isfile(IMAGE) && isfile(STAMP) || return nothing
    julia, checksum = split(read(STAMP, String), '\n')
    return isfile(julia) && checksum == image_checksum() ? String(julia) : nothing
end
