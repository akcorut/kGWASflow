# Changelog

## [1.2.0](https://github.com/akcorut/kGWASflow/compare/v1.0.0...v1.2.0) (2023-04-17)


### Features

* add a function to extract chromosome names and lengths from a SAM file ([413ead9](https://github.com/akcorut/kGWASflow/commit/413ead9fdf373c4809c2ade506b5d77a82821634))
* add a rule for genome indexing ([f2d596e](https://github.com/akcorut/kGWASflow/commit/f2d596e8d21af1be92a78d68b1a2ff6d7efb6eb3))
* add a rule to convert BAM files to BED ([a948d98](https://github.com/akcorut/kGWASflow/commit/a948d987de60dbbc7615cfbb7618a0192217d9e5))
* add a rule to convert kmers table to PLINK ([8c4bfad](https://github.com/akcorut/kGWASflow/commit/8c4bfad582c84d510aef8a81be20536a572a89df))
* add error handling ([da92bdb](https://github.com/akcorut/kGWASflow/commit/da92bdbb749ab7b841cd87fe6ba68f45627d2a3e))
* add k-mer counts dist. plots ([7df4787](https://github.com/akcorut/kGWASflow/commit/7df4787221c8f5d25cd0ee8e6c8689a25ae80686))
* add k-mer counts dist. plots ([79427b4](https://github.com/akcorut/kGWASflow/commit/79427b447dbba2b6ea20ccd519574785f509b351))
* add k-mer counts dist. plots ([c45c412](https://github.com/akcorut/kGWASflow/commit/c45c4125101a1526952877bee146254c17ca341c))
* add kGWASflow.py ([698dcb8](https://github.com/akcorut/kGWASflow/commit/698dcb83ecfa4d160b793adf4e844216b98db7f5))
* add new functions, update checkpoints + target outputs ([6dac404](https://github.com/akcorut/kGWASflow/commit/6dac404e0c2f94928203ae03bd688a851537ddba))
* add rules for converting BAM to BED ([33ed058](https://github.com/akcorut/kGWASflow/commit/33ed058e90d017f0cd60932f65494237c848eeae))
* add rules for generating igv reports (igv-report) ([2c4b4fe](https://github.com/akcorut/kGWASflow/commit/2c4b4fe79cf78a3e3ed96b6fa3d8f84509079df0))
* improve the kmers GWAS summary reports ([9d5de9e](https://github.com/akcorut/kGWASflow/commit/9d5de9e60e279c585cbd30d2cdec4d157d6ccd6c))
* udpate config files + add new parameters ([38c939b](https://github.com/akcorut/kGWASflow/commit/38c939ba9cac5d8f3a337a5366116f02c7111d61))
* update blast wrapper to v1.25.0 ([3da03ba](https://github.com/akcorut/kGWASflow/commit/3da03bad243d15c07dde8f9a15c2c8b490900fdd))
* update to kmc 3.2.1 ([c933fda](https://github.com/akcorut/kGWASflow/commit/c933fdafb5ba16c42d9becd15464960280485cb5))


### Bug Fixes

* change output, log paths and report caption ([827dcad](https://github.com/akcorut/kGWASflow/commit/827dcadb28cd81cb35853449c8e07f0913cd7bdb))
* change x and y-axis ticks and labels ([c7cdf87](https://github.com/akcorut/kGWASflow/commit/c7cdf87ff1246b8cc3374a8785afdeab987b7e4c))
* fixed typo ([1238463](https://github.com/akcorut/kGWASflow/commit/12384631520ceba8ca6a9c37d7995a884994783a))
* make the script compatible to new kmers GWAS summary reports ([3dbecf0](https://github.com/akcorut/kGWASflow/commit/3dbecf0f5661b4c641e37476dbdbe587f1bbe710))
* specify the spades version ([35d5690](https://github.com/akcorut/kGWASflow/commit/35d5690f113222c168907f4ef2623bcac5d256d4))
* update input file paths ([05d7e6b](https://github.com/akcorut/kGWASflow/commit/05d7e6b8499920b3c175d7def7909dc241b27adc))
* update input, ouput and params paths ([4637ed6](https://github.com/akcorut/kGWASflow/commit/4637ed6349755397a966bba7e5b4666b03ecd99b))
* update output directory paths + add touch ([a0500c2](https://github.com/akcorut/kGWASflow/commit/a0500c26248408f288acaa984f80447b5020c166))


### Performance Improvements

* update snakemake version to 7.25.0 ([5b999c4](https://github.com/akcorut/kGWASflow/commit/5b999c48382ddbf647ef235ee214783fee572bcb))
* update to fasterq-dump wrapper v1.23.5 ([373b406](https://github.com/akcorut/kGWASflow/commit/373b4063c2c13d1485de77863b0bb6f5f6296492))
