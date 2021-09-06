# Stochastic Generation of (t,s) Sample Sequences

This is a C++ implementation of a large portion of the techniques and algorithms described in ["Stochastic Generation of (t,s) Sample Sequences"](https://www.seanet.com/~myandper/abstract/egsr21.htm), by Helmer (that's me!), Christensen, and Kensler (2021). 

Short-version for people who aren't familiar with the paper: if you want to generate Owen-scrambled Sobol', Halton, or Faure sequences, or pmj02(bn) sequences, I don't think there are any faster or simpler open-source implementations out there.

## Generating Samples

This repository has two main goals. The first is to make it easy for people to generate various sequences from the paper, with different options that the techniques provide. The second is to provide reference implementations, so that others can build on these sequences and techniques and maybe better understand the paper.

To generate samples, one should first build the "generate_samples" utility. The example commands given here will apply to shells in Linux or macOS. In a shell, you'd first run:

<pre><code>make generate_samples</code></pre>

One can then generate different sequences by using the --seq= flag to specify what kind of sample sequence to generate, --n= to specify the number samples, and --nd= to specify the number of dimensions of a multidimensional sequence. For example, one can generate 4096 samples of a 5D stochastic Sobol' sequence with the following Bash command:

<pre><code>./generate_samples --seq=ssobol --n=4096 --nd=5 > /my/samples/directory/ssobol.txt</code></pre>

### Available Sample Sequences

This library implements seven sample sequences (or algorithms) that can be generated:

#### ssobol

The stochastic Sobol' sequence. One can generate up to 64 dimensions. The first two dimensions form a base-2 (0,2)-sequence, which means that they are stratified on all base-2 elementary intervals for a power-of-two prefix of samples. This is shown for the first 16 samples:

<p align='center'>
	<img src='https://github.com/Andrew-Helmer/stochastic-generation/blob/main/images/ssobol02.jpg'><br>
</p>

Actually this is true for all power-of-two "disjoint subsequences". So it's not just the first 16 samples, but the next 16, and the 16 after that, and so on.

#### pmj02

We can generate the pmj02 sequence from "Progressive Multi-Jittered Sample Sequences" by Christensen, Kensler, and Kilpatrick (2018), but using the Stochastic Generation approach as described in our supplemental material. This is only a 2D sequence, so the flag --nd=2 must be set. For most intents and purposes, this will be the same as the first two dimensions of the ssobol sequence, although the points will be ordered differently within each power of two.

#### sfaure(03|05|07|011)

This library implements sfaure03, sfaure05, sfaure07, and sfaure011 sequences. The sfaure03 sequence is a base-3 (0,3)-sequence. This means that it's a 3-dimensional sequence where any power-of-three disjoint subsequence of samples is going to be distributed on all base-3 elementary intervals. Here's Figure 7 from the paper, showing how well stratified that is, in 3D, in all three 2D projections, and all three 1D projections:

<p align='center'>
	<img src='https://github.com/Andrew-Helmer/stochastic-generation/blob/main/images/sfaure03.jpg'><br>
</p>

If you wanted to estimate a 3D integral, and you can accept the limitation of only using powers of three, this will give very low error. The same is true using an sfaure05 sequence for a 5D integral, at powers of five. The code would be very easy to extend to higher prime bases (e.g. a stochastic Faure (0,13)-sequence), but they probably wouldn't be very useful in practice. 

#### shalton

The stochastic Halton sequence. The Halton sequence doesn't have all of the elementary interval stratifications of the sfaure sequences or the first two dimensions of the ssobol sequence, but it does a pretty nice job of guaranteeing stratifications on all lower-dimensional projections across different sample counts. 

The Halton sequence is constructed with subsequent prime number bases for each dimension, e.g. the first two dimension is base-2, the second dimension is base-3, the third is base-5, base-7, etc. This implementation goes up to 10 dimensions (base 29), though it could be easily extended to higher.

So for any subsequence (not just disjoint!) of 2<sup>a</sup> * 3<sup>b</sup> * 5<sup>c</sup> ... points, where the exponents are >= 0, the points are stratified in the 2<sup>a</sup> x 3<sup>b</sup> x 5<sup>c</sup> ... grid within the unit hypercube. Those grids include all lower dimensional projections (i.e. when some exponents are zero). For example, any subsequence of 72 points will be stratified in the 8x9 grid in the first two dimensions. Any subsequence of 15 points will be stratified in the 3x5 grid in the second and third dimension. And so on.

### Best-Candidate Sampling

Every sequence can be augmented with best-candidate sampling in the first two dimensions using the "--bn2d" flag. Rather than picking any random point in a valid strata, many possible candidates are evaluated, and the one that is furthest from all previous points (in those two dimensions) will be used. This can slightly improve point spacing on these sequences. For instance, to generate a pmj02bn sequence, the command would be:

<pre><code>./generate_samples --seq=pmj02 --n=4096 --nd=2 > /my/samples/directory/pmj02bn.txt</code></pre>

Here's a simple comparison of 64 pmj02 points vs. 64 pmj02bn points. The spacing is slightly improved overall, though not dramatically.

<p align='center'>
	<img src='https://github.com/Andrew-Helmer/stochastic-generation/blob/main/images/pmj02bn.jpg'><br>
</p>

This was limited to only the first two dimensions for simplicity. But it would be easy to extend the implementation to apply to other dimensions, where it would be appropriate. It could be useful for higher dimension pairs of the Faure or Halton sequences. If you wanted to do this, take a look at the function "GetBestFaurePoint" in sfaure.cpp, and "GetBestHaltonPoint" in shalton.cpp.

### Correlated Swapping vs. Owen-Scrambling

By default, the shalton and sfaure sequences are generated with correlated swapping, which is described in Section 4.7 of the paper. You can instead use "uncorrelated swapping", which would generate sequences that are identical to fully Owen-scrambled sequences, using the --owen flag. For example, the command:

<pre><code>./generate_samples --seq=sfaure05 --owen --n=3125 --nd=5 > /my/samples/directory/sfaure05.txt</code></pre>

will generate an Owen-scrambled Faure (0,5)-sequence.

We haven't observed much differences between the quality of sequences generated with correlated swapping and Owen-scrambling, but it's still nice to include both here and provide reference implementations for either. One thing to note is that the choice of strata is done randomly even when best-candidate sampling is used. One can get a better sequence when also incorporating multiple strata choices into the best-candidate sampling, and we did this for Figure 2c of the supplemental material. This was omitted because the code to keep track of the available strata is quite complex and ugly, and readability was an important goal for this code.

### Shuffling and Decorrelation

With the --shuffle flag, the ssobol and sfaure sequences will be *progressively* shuffled, which will decorrelate subsequent runs of a particular sequence. Because the shuffling is done progressively, any power-of-b prefix of the sequence will be unchanged, so best-candidate sampling will still improve distances. See Section 5.2 of the paper for more information on this.

## Performance

This repository is meant to somewhat balance generality and readability with the high-performance aspects of the paper. This means that the performance is not as fast as it would be, if the code was written for a single specific sequence and set of options. For example, Listing 2 from the paper, which shows a simple implementation of *only* a 2D Sobol' (0,2)-sequence, without best-candidate sampling, shuffling, or the possibility of higher dimensions, will be about 3x faster than generating the same sequence using this code. Even so, the code in this library is faster than any other published method.

We also didn't include the "index mapping" precomputation of the sfaure sequences. For one-off generation of each sequence, it wouldn't be any faster, so it didn't make a ton of sense to include here. However if you wanted to generate a lot of sfaure sequences, that precomputation would be helpful.

## Progressive Shuffling Indices

Also included is a command-line utility to generate shuffled index arrays. You would *make* this separately:

<pre><code>make shuffle_indices</code></pre>

And then if you wanted to generate an array of base-3 shuffled indices, of length 9, you would do something like:

<pre><code>./shuffle_indices --n=9 --b=3</code></pre>

which might output:

<pre><code>4 5 3 0 2 1 7 8 6</code></pre>

If you set --progressive, this will do progressive shuffling, maintaining the original prefix samples at each power of the base. For example

<pre><code>./shuffle_indices --n=9 --b=3 --progressive</code></pre>

might give:

<pre><code>0 1 2 3 5 4 8 7 6</code></pre>

## Xor-value Generation

The xor_values.py file can be used to generate the xor-values for the sfaure and pmj02 sequences, straight from their mathematical definitions. The xor-values for the Sobol' sequence are based on the choice of direction numbers, which is itself a complex process, so this was omitted. But the PBRT v3 renderer has a set of generator matrices for the Sobol' sequence, which were used to obtain the xor-values in this code. You can find PBRT's generator matrices here: https://raw.githubusercontent.com/mmp/pbrt-v3/main/src/core/sobolmatrices.cpp

## Licensing

Licensed under the MIT Open Source license. See the [LICENSE](/LICENSE) file. 
