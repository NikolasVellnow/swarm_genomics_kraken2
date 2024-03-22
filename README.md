# Analysis of unmapped reads with kraken2 #


There is often interesting biological information ...blalba

We assume we have miniconda installed, and we have a miniconda environment called "kraken" in which kraken2 is installed.


First we activate the miniconda environment (named "kraken") in which we have installed the kraken software so we can use it:

```sh
conda activate kraken
```

Now we first have to download and build a database containing the *l*-mers and the taxonomic information. We are using the *kraken2-build --standard* command for that. This command will download and build the standard Kraken 2 database. Afterwards we can use that database to classify sequencing reads as belonging to specific taxa. You need to replace $DBPATH with the path to where you want to save the database and the name you want to give the database. You can also choose the number of threads you  want kraken2 to use. Just replace $THREADNUM with the number of threads you want to use, e.g. 20. The database building will take a while.

```sh
kraken2-build --standard --threads $THREADNUM --db $DBPATH
```

But we can also build a custom database ourselves. For this, we first download the taxonomy information and the information for the different taxa we want to have in our database:

```sh

# download taxonomy
kraken2-build --download-taxonomy --threads $THREADNUM --db $DBPATH

# download most dbs
echo 'archaea'
kraken2-build --download-library archaea --threads $THREADNUM --db $DBPATH

echo 'bacteria'
kraken2-build --download-library bacteria --threads $THREADNUM --db $DBPATH

echo 'plasmid'
kraken2-build --download-library plasmid --threads $THREADNUM --db $DBPATH

echo 'viral'
kraken2-build --download-library viral --threads $THREADNUM --db $DBPATH

echo 'human'
kraken2-build --download-library human --threads $THREADNUM --db $DBPATH

echo 'fungi'
kraken2-build --download-library fungi --threads $THREADNUM --db $DBPATH

echo 'plant'
kraken2-build --download-library plant --threads $THREADNUM --db $DBPATH

echo 'protozoa'
kraken2-build --download-library protozoa --threads $THREADNUM --db $DBPATH

echo 'UniVec_Core'
kraken2-build --download-library UniVec_Core --threads $THREADNUM --db $DBPATH


```

Then we build the database:


```sh
kraken2-build --build --db $DBPATH --threads $THREADNUM

```


However, the easiest way is to just download an already made database from the this page: https://benlangmead.github.io/aws-indexes/k2

For example, we could download the PlusPFP-16 database (just click on the ".tar.gz" in that row). After downloading it we need to unzip the tar ball:

```sh

tar -xvzf k2_pluspfp_16gb_20240112.tar.gz

```

We can now classify the reads from the file of interest using Kraken 2. Let's say we want to classify the unmapped reads from whole-genome sequencing of a great tit.


```sh





kraken2 \
--db ${DBPATH} \
--output ${OUTPUTNAME} \
--use-names \
--report ${REPORTNAME} \
--classified-out ${CLASSIFIEDNAME} \
--unclassified-out ${UNCLASSIFIEDNAME} \
--confidence 0.1 \
--threads ${THREADNUM} \
${INPUTNAME}

```








