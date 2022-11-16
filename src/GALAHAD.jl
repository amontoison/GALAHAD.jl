module GALAHAD

galahad_version() = v"5.0.0"

# packages without dependencies.
include("C/bsc.jl")
include("C/convert.jl")
include("C/fit.jl")
include("C/glrt.jl")
include("C/gls.jl")
include("C/gltr.jl")
include("C/hash.jl")
include("C/hsl.jl")
include("C/ir.jl")
include("C/l2rt.jl")
include("C/lhs.jl")
include("C/lms.jl")
include("C/lsrt.jl")
include("C/lstr.jl")
include("C/presolve.jl")
include("C/roots.jl")
include("C/rpd.jl")
include("C/scu.jl")
include("C/sec.jl")
include("C/sha.jl")
include("C/sils.jl")
include("C/ugo.jl")
include("C/ssids.jl")

# sls requires sils, ma57, ma77, ma86, ma87, ma97, ssids, mc64, mc68.
include("C/sls.jl")

# rqs requires sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, ir.
include("C/rqs.jl")

# dps requires sls, sils, ma57, ma77, ma86, ma87, ma97, ssids, mc64, mc68.
include("C/dps.jl")

# psls requires sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, mi28.
include("C/psls.jl")

# arc requires rqs, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, ir, glrt, dps, psls, mi28, lms, sha.
include("C/arc.jl")

# trs requires sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, ir.
include("C/trs.jl")

# trb requires trs, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, ir, gltr, psls, mi28, lms, sha.
include("C/trb.jl")

# bgo requires trb, trs, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, ir, gltr, psls, mi28, lms, sha, ugo, lhs.
include("C/bgo.jl")

# uls requires gls, ma48.
include("C/uls.jl")

# sbls requires sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48.
include("C/sbls.jl")

# blls requires sbls, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, convert.
include("C/blls.jl")

# bqp requires sbls, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48.
include("C/bqp.jl")

# fdc requires sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48.
include("C/fdc.jl")

# cro requires sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, sbls, uls, gls, ma48, ir, scu.
include("C/cro.jl")

# bqpb requires fdc, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, sbls, fit, roots, cro, ir, scu, rpd.
include("C/bqpb.jl")

# ccqp requires fdc, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, sbls, fit, roots, cro, ir, scu, rpd.
include("C/ccqp.jl")

# clls requires fdc, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, sbls, fit, roots, cro, ir, scu, rpd.
include("C/clls.jl")

# cqp requires fdc, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, sbls, fit, roots, cro, ir, scu, rpd.
include("C/cqp.jl")

# dgo requires trb, trs, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, ir, gltr, psls, mi28, lms, sha, ugo, hash.
include("C/dgo.jl")

# dqp requires fdc, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, sbls, gltr, scu, rpd.
include("C/dqp.jl")

# eqp requires fdc, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, sbls, gltr.
include("C/eqp.jl")

# lpa requires rpd.
include("C/lpa.jl")

# lpb requires fdc, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, sbls, fit, roots, cro, ir, scu, rpd.
include("C/lpb.jl")

# lsqp requires fdc, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, sbls.
include("C/lsqp.jl")

# nls requires rqs, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, ir, glrt, psls, mi28, bsc, roots.
include("C/nls.jl")

# qpa requires sls, sils, ma57, ma77, ma86, ma87, ma97, ssids, mc64, mc68.
include("C/qpa.jl")

# qpb requires lsqp, fdc, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, sbls, gltr, fit.
include("C/qpb.jl")

# slls requires sbls, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, convert.
include("C/slls.jl")

# tru requires trs, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, ir, gltr, dps, psls, mi28, lms, sec, sha.
include("C/tru.jl")

# wcp requires fdc, sls, sils, ma57, ma77, ma86, ma87, ma97, ssids,
# mc64, mc68, uls, gls, ma48, sbls.
include("C/wcp.jl")

end # module GALAHAD
