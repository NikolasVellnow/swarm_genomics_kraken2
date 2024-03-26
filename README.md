# Analysis of unmapped reads with kraken2 #


There is often interesting biological information left to explore in the reads that did not map to a reference genome (e.g., Chrisman et al. 2022). This document will show you the first steps of how to do such an anylsis of unmapped reads using the software Kraken 2.

Here, we assume that we have miniconda installed, and that we have a miniconda environment called "kraken" in which Kraken 2 and KrakenTools are installed.

First we activate the miniconda environment (named "kraken") in which we have installed the kraken software so we can use it:

```sh
conda activate kraken
```

## Building a database ##
Now we first have to download and build a database containing the *l*-mers and the taxonomic information. We could (i) download an already built database from https://benlangmead.github.io/aws-indexes/k2, (ii) use the *kraken2-build --standard* command, or (iii) build our own custom database with all the taxa we want.
For an easy start into the usage of Kraken 2 I will describe how to use a pre-build database, i.e. (iii). Note, however, that it might be necessary to build your own custom database to get the best and most meaningful results.

First, we just download an already made database from the this page: https://benlangmead.github.io/aws-indexes/k2
For example, we can download the PlusPFP-16 database (just click on the ".tar.gz" in that row). After downloading it we need to unzip the tar ball:

```sh
tar -xvzf k2_pluspfp_16gb_20240112.tar.gz
```
Now we have the database we can compare our unmapped reads to.

## Classifying reads ##
As an example, I will now classify the unmapped reads from whole-genome sequencing project of a great tit (*Parus major*). Form this we use the *kraken2* command.
You need to replace $DBPATH with the path to where you have saved the database and the name you gave the database. You can also choose the number of threads you  want kraken2 to use. Just replace $THREADNUM with the number of threads you want to use, e.g. 20. The flags *--classified-out* and
*--unclassified-out* are optional.

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

This command will take a while and it will give you a report and a kraken-output. The output (.kraken ending) lists all the reads and as what they were classified. The Report (.kreport ending) is the most useful one, because you can see how many reads were classified as what taxon.

## Visualizing the results ##
We can now use the script *kreport2krona.py* from KrakenTools to make a nice visualization of the results.

```sh
python3 kreport2krona.py -r REPORTNAME -o SAMPLE.krona
ktImportText SAMPLE.krona -o SAMPLE.krona.html
```
## Resources ##
You can learn more about the use of Kraken 2, e.g. how to build your own custom database, here: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
And here is information about KrakenTools: https://github.com/jenniferlu717/KrakenTools/blob/master/README.md



