## dBGMapper  （building version）

Author

```
~Changyong Yu (Northeastern University in CHINA)
~Chu Zhao (Northeastern University in CHINA)
```

### 1.Introduction

dBGMapper is a tool building for bioinformatics read mapping, which based on the seed-and-extension strategy. It works on the dBG(deBruijn Graph) data generated from [StLiter](https://github.com/BioLab-cz/StLiter) and designed to find all matches under the preset threshold. 

Using dBG as index data can solve high repetitiveness in gene sequence and therefore has excellent performance in the use of space. We utilize pigeonhole as seed division principle and extend the candidate position of seed. FM-index is used for the very beginning of finding the seed exact match position in reference genome. Then use kmer information appeared in dBG as basic extension and alignment data. Taking various extension strategies according to the frequency of seed candidates. Finally, all candidate positions are calculated and used to select the correct or the best match of read mapping.

---

### 2.Schedule

1. ~~dBG attribute analyses with different kmer lengths for memory usage feasibility analysis.(**Accomplished**)~~

2. ~~dBG data input and appropriate data structure organization for quick sort and quick query.(**Accomplished**)~~

3. ~~obtain of frequency and precise location to short gene sequence.(**Accomplished**)~~

4. ~~location in dBG data with given kmer information.(**Accomplished**)~~

5. ~~kmer extension in the dBG according to above-mentioned location.(**Accomplished**)~~

6. early termination in extension strategy to optimize above-mentioned process.(**building**)

   

### 3.To be continued

