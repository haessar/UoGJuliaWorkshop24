# The (Newton)-Raphson-Simpson method

This exercise looks at simple function construction and plotting in Julia, with some use of the type system and lists.

## Introduction

$x_{n+1} = x_n - \frac{f(x_n)}{f′(x_n)}$

Commonly known as "Newton's method", despite the most significant developments in its modern form being due to Joseph Raphson and Thomas Simpson (and Newton's version being independently developed in historically and geographically disparate locations, from the Middle East to Japan), this iterative method of finding the roots of a function `f`, given its known derivative `f'` should be familiar to most of you.

In this exercise, we will explore the implementation of the (Newton)-Raphson-Simpson method - henceforth "NRS" - in Julia, as a way to experience many of the features of the language.

## Getting Started

If you haven't already, launch a `julia` REPL - or a Jupyter notebook with a Julia kernel. 

For the first parts of this exercise, you will not need to load any additional packages, everything is just core Julia.

Firstly, we'll write a simple function to represent the delta between one estimate for the root, $x_n$ and the next, $x_{n+1}$, given the function and its derivative. 
Remember, in Julia you can type most special characters by typing their $\LaTeX$ representation, and then pressing tab.

So, our "nrs_delta" can be defined as, typing into the REPL, and pressing enter:

```julia
"""nrsδ(x,f,f′)
calculate the next correction to x, our estimate of a root of f, with derivative f′
"""
nrsδ(x, f, f′) = - f(x) / f′(x)
```

where you can type $\delta$ with `\delta` and then `tab`, and the prime on the f' as `\prime` and then `tab` . Here the `""" """` quoted string at the start provides the help text for anyone querying the function - in the REPL, you can use `<ALT>+<ENTER>` (in Linux), `<Opt>+<ENTER>` (MacOS), or `<ALT>+<SHIFT>+<ENTER>` at the end of the string to make an "escaped newline" so you can continue writing input.

* Test that the help text works by using the help mode (with ?) to query help for our function's name.

We'll also define some test `F` and `dF` functions so we can easily test our code [^1]: 

```julia
F(x) = x^3 - 8
dF(x) = 3x^2
```

notice, also, how Julia lets us write multiplication by juxtaposing a number and a variable, if it isn't ambiguous.

### NRS function, on Numbers

The NRS algorithm itself needs to repeatedly evaluate new deltas on our estimates, and exit when the magnitude of the delta is small enough - say, 1e-6 

To implement this as a function, it needs to take three things: an initial $x_0$ - ideally some kind of `Number` - and the function and its derivative.

Write a function called `nrs` which implements the following logic:

* take in three things:
    * $x_0$, which must be a kind of `Number`
    * $f$
    * $f^\prime$

($f$ and $f^\prime$ are both functions, but don't try to specify that in the function signature!)

Also write a suitable docstring for the function. 

Function body:
* $\epsilon \leftarrow 1\times10^{-6}$, our closeness criterion
* $x_n \leftarrow x_0$ 
* set a local variable $\Delta$ equal to 1.0 (so we have a placeholder value

* loop, if $\Delta > \epsilon$
    * $\Delta \leftarrow nrs\delta(x_n,f,f^\prime) $
    * $x_n \leftarrow x_n + \Delta$
* return the value of $x_n$ when we exit the loop

You can use the `abs` function to get the absolute value of a number.

* Try this new function out with various values - it should work with any kind of `Number`, from integers through to even complex numbers.


### Extension - Limits

At present, this function is flawed - if a particular value for $x_0$ doesn't lead to a converging solution, it will iterate for ever.

By adding a new variable inside our function, called `count`, try modifying the NRS function we just wrote to also stop if it loops more than 40 times.

(You can either add a second condition to the `while` loop, or use `break` inside the loop to break out).

If you are in the REPL, you can add new lines to a multi-line entry in the history by pressing `ESC` before pressing `enter`. (Just pressing `enter` will execute the multi-line entry instead).

## Extending the NRS function to store history

Currently, our NRS algorithm is perfectly fine, but it might be nice to have an alternative version which takes a Vector (a 1-d Array) with an initial element representing $x_0$, and appends the history of all $x_n$ to it.

That way, we have a record of the path the function has taken to converge to the root it finds.

We'll define a new method for the NRS function that does this - multiple dispatch will select the right version for us based on the type of the first argument.

We'll need to make only a few changes relative to the first method for NRS:

**Firstly** - this method's first argument is a Vector of values of type T, where T must be some kind of Number. We can express this using the "where" clause in our function definition, after the parameter list.
e.g
```julia
function myfunc(v::Vector{T}) where T <: Number
```
We should call our new first parameter `xs`, where its first element will be effectively $x_0$

**Secondly**, rather than `x_n` being assigned to from `xs` directly, it needs to take the first element of `xs`. 

**Thirdly**, we need to add each new value of `x_n` to the `xs` vector. Using help mode (`?`) investigate the `push!` method, which can be used to implement this.

**Finally**, we need to return `xs` and not `x_n` at the end of the function.

### Testing it out

* Make a version of NRS with the above changes. 

Test it out: try passing the value `5.0` to NRS, and check it still just returns a single value. Then try passing `[5.0]` and see what you get as a result!

### An Aside

*Technically*, this method will actually modify the Vector we pass to it, as Julia implements *pass by sharing*. 

The convention in Julia is that functions that modify their arguments - or other state - should have a `!` at the end of their name. So, this should really be a new function called `nrs!`. 

* Make a version of `nrs` which does not modify its Vector argument - making a deep copy of it instead, and returning that.


## C) Plotting!

This final part of the exercise needs the `Plots` library.

If you haven't installed this before, then firstly we'll need to go in to `Pkg` mode to set it up.

### Pkg in the REPL

If you're in the REPL, press `]` to go into Pkg mode, and then type

```Pkg
    add Plots
```

and wait for Pkg to fetch and precompile all the prerequisites for the Plots package for you.

### Pkg in Jupyter notebooks

Pkg mode doesn't work in notebooks, so instead you will need to load the Pkg library, and use its API.

type

```julia
    using Pkg
    Pkg.add("Plots")
```

and again wait for Pkg to sort out your dependencies.

### Plots

Now we've got `Plots` installed, we need to load it with

```julia
    using Plots
```

(there will be a brief delay at this point as Plots sets itself up)

If you're in a Jupyter notebook, you will get "inline" plots after each cell. If you're in a REPL on a machine with graphical display, Plots will instead open you a separate window to display figures.

#### Line plots

Let's start by displaying the function that we're finding zeros for.  We'll stick to the Real line for now, so we can use a simple `line` plot, which is also the default.

We know, of course, that the Real zero for our function is at $x = 2$, so lets set up an x range from 0 to 3 to get a view of it.

`Plots` can take lots of arguments to set up its x and y ordinates. For us it's easier to specify a range for x, and a function for y (`Plots` will call the function on all the x values it needs to)

A range in Julia looks like `start:step:end`, so...

```julia
    plot(0.0:0.01:3.0, F)
```

You should get a nice, unsurprising, plot of our function $x^3 - 8$ from 0 to 3.

#### Scatter plots

Now lets modify our figure to add in a trace of all the points we evaluate with our NRS function.

We can collect our history of points with our enhanced version of the nrs function (for a given $x_0$, here $1.2$ )

```julia
    history = nrs([1.2], F, dF)
```

Now we have a bunch of x values, but we also need the values of F(x) for each of them. We can either broadcast over the history vector, or simply pass our plotting method the function to use to generate the y values directly.

* using one of the mutating plot methods, add a scatter plot to the figure showing the points on the curve that we evaluate.

*Remember* you can use the ? help mode to find out how to use functions. You can probably guess the name of the function that plots scatter plots, and that the mutating version ends with a particular character.

### Going further: the complex plane

This is nice so far, but we know fine well that this problem is well-defined over the entire complex plane. 

In order to display a reasonable representation of a 4 dimensional space (the 2 components of our complex-valued function at each point in the complex plane), we'll stick to simply plotting the absolute value of F at any point.

`contour` is a good function for this representation - it generally wants 3 values for x,y,z, but you can use the help mode to check how it works.

x and y can just be ranges as before (lets go from -3 to 3 in both directions to get a good view).

z will need to be a function that maps (x,y) back to complex values, calls F and then calls abs on that. We should write this as an anonymous function.

* make the appropriate contour plot of F 

If we want to display the trace of an NRS in the complex plane, we can just forgo representing the value of F(x) and plot the x,y coordinates of each $x_n$ - using `real()` and `imag()` as broadcasted functions to get the x and y components of our points.

* collect a history of points starting with an *complex* value for $x_0$

* overlay an appropriate scatter plot on the contour.

## Further work

At present our scatter plot makes it hard to see the ordering of the iteration. 

1. We could provide an array to the `markercolor` property of our scatter plot, allowing each point to take a different colour.

1. We've also only really tested this with one F and dF - try passing different candidate functions to NRS and see how different initial values converge!

### Advanced extension

These functions currently all require the user to provide both `f` and `f′`, which means that they are prone to user error.

There are many automatic differentiation packages available for Julia - such as `Zygote` - listed at [Julia Diff](https://juliadiff.org), which we could use to find `f′` directly and efficiently.

Install `Zygote`.

Using its `gradient` method, and the fact that (for a holomorphic function), the complex gradient at $z$ is

${(\frac{d\Re(f)}{dz}\bigr\rvert_z) }^*$
where `conj` is the complex conjugation operator in Julia

* make versions of `nrsδ` and `nrs` that take only the initial value and the function to find the root of (assuming it is holomorphic).


[^1]: Technically, these functions may not be as fast as they could be, as the constant 8 will not always have the same type as the value x passed to F (or dF).
We can write a more complex version of the same function which could be faster as:
```julia
F(x::T) = x^3 - T(8)
```
where here we can use the name of the type of x, T, as a function to coerce the value 8 to the same type explicitly.