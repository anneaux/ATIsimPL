
## Brief explanation of the arguments of the `get_thimble()` function

```

get_thimble(S::Function, drv::Function, tmin::Float64, tmax::Float64;
    Nflow::Int64 = 60,
    Δinit::Float64 = 10.,
    flowstepfactor::Float64 = 2.,
    h_threshold::Float64 = -300.,
    gradnthreshold::Float64 = 1.,
    subdividethreshold::Float64 = 4.
    )

``` 

### in short
- `Nflow`: how many flow steps to do
- `Δinit`: how fine to chop the integration domain in the beginning
- `flowstepfactor`: how big the flow steps are
- `h_threshold0`: how deep down into the valley the contour goes
- `gradnthreshold`: how precise we are around the saddle points (ignore)
- `subdividethreshold`: how fine we chop the integration domain into pieces (smaller number --> finer)

### slightly longer story

The integration domain (here a line along the real time axis) is represented as a set of points and the resulting line segments (with distance < `Δinit`).
In each iteration step of the flow procedure, each of the points "flow" in the direction of the normalised gradient* (see step in the `flow()` function [here](https://github.com/anneaux/ATIsimPL/blob/5e00e120d85b188607c55683feba15b7c4e9fdb1/PLIntegration1D.jl#L95C1-L95C67 ), multiplied with the `flowstepfactor`.
If the distance between two points is too big (greater than `subdividethreshold`), we insert a new one in the middle.
If the imaginary part of the phase function drops below a certain threshold (namely, `h_threshold`), the point is turned inactive and won't flow any longer. That is the part of the contour flows deep down into the valley where the integrand vanishes.

The algorithm is adapted (and nicely explained) on [Job Feldbrugge's website](https://p-lpi.github.io/).

(*) the gradient will be normalised if its norm is greater than `gradnthreshold`.



