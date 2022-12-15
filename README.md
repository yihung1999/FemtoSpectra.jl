# FemtoSpectra.jl

## How to Use

```julia
using FemtoSpectra
```

Create input & output power data `(iDat, oDat)`:

```julia
iDat = [472.0, 552.0, 692.0, 909.0, 1146.0, 1483.0]
oDat = [100.0, 122.5, 147.2, 193.1, 243.2, 320.8]
```

Then, create the converter (based on linear regression):

```julia
powerIO = PowerIO(; xdat=iDat, ydat=oDat)
```

Check two parameter (slope, intercept):

```Julia
powerIO.wVec
```

Enter numbers in brackets
EX: 114.4
```Julia
powerIO(114.4)
```