## dBGMapper  （building version）

Author

```
~Changyong Yu (Northeastern University in CHINA)
~Chu Zhao (Northeastern University in CHINA)
```

### 1.Introduction

dBGMapper is a tool building for bioinformatics read mapping, which based on the seed-and-extension strategy. It works on the dBG(deBruijn Graph) data generated from [StLiter](https://github.com/BioLab-cz/StLiter) and designed to find all matches under the preset threshold. 

Using dBG as index data can solve high repetitiveness in gene sequence and therefore has excellent performance in the use of space. We utilize pigeonhole as seed division principle and extend the candidate position of seed. FM-index is used for the very beginning of finding the seed exact match position in reference genome. Then use kmer information appeared in dBG as basic extension and alignment data.

---

### 2.Schedule

