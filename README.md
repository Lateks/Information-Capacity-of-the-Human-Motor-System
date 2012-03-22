# Scripts for information capacity evaluation

## Information capacity evaluation scripts
- `eval.R` (the actual functions to call when you want to evaluate data)
- `infocapacity.R` (residual calculators, the actual evaluation functions, normalization etc.)
- `data_handling.R` (handling of files, aligning data etc.)

### The actual functions to call:

The function `residual_complexity(fps, pca)` calculates throughputs for the whole directory, assuming the original files (named `01.txt`, `02.txt` etc.) are in the current directory and a subdirectory `alignment` contains the alignment vectors indicating frame duplications as produced by CTW (named `1_ali_2.txt` etc.).

The function `pair_residual_complexity(seqnum_a, seqnum_b, fps, pca)` calculates throughputs for the given pair of sequences, indicated by the given sequence numbers.

The function `subdir_based_residual_complexity(fps, pca)` calculates throughputs for the current directory assuming that each subdirectory `01`, `02` etc. contains all the repetitions of one type of sequence. The directory structure and file naming is assumed to be like this:

    - working directory
      |
      -- 01 -- contains three repetitions of sequence #1
      |  |
      |  * 01.txt
      |  * 02.txt
      |  * 03.txt
      |  -- alignment -- contains the alignment vectors
      |     |
      |     * 1_ali_2.txt
      |     * 1_ali_3.txt
      |     * ...
      -- 02 -- contains two repetitions of sequence #2
      |  |
      |  * 01.txt
      |  * 02.txt
      |  -- alignment
      |     |
      |     * 1_ali_2.txt
     ...    * 2_ali_1.txt

Of the parameters for each of these functions, `fps` indicates the framerate in frames per second (default is 120) and `pca` is a logical value indicating that PCA should be used (default is false, PCA is not performed).

## Plotting scripts
- `plot_aligned.R` (plotting aligned sequences into pdf files)
- `plot_outliers.R` (diagnostic functions for "outlier" residuals, very messy prototypes)
