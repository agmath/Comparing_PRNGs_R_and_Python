# Comparing p-RNGs in R and Python

A repository for code to accompany the paper titled A Comparison of pseudo-Random Number Generators in R and Python by Mata-Toledo, Dumnich, Ryan, and Gilbert. Functions utilized and referenced throughout the paper are contained in the `CurrentPaperCode` folder. Any script names beginning with `generate_*` are helper functions referenced by the other scripts in the folder. These helper functions are used to generate the random sequences of values.

In our paper, we explored the implementations of the Mersenne Twister algorithm in both R and Python. We also considered the Permuted Congruential Generator family and the XORshift128+ algorithm. We've organized a `results.qmd` file to allow interested readers to reproduce all of the results from our paper. You can find that file in the `CurrentPaperCode` directory or by following [this direct link to the file](https://github.com/agmath/Comparing_PRNGs_R_and_Python/tree/main/CurrentPaperCode/Results.qmd). Running this file depends on functionality also contained in the `CurrentPaperCode` directory. Cloning the repository is the quickest method to begin reproducing our results.

