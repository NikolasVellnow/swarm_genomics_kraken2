# Analysis of unmapped reads with kraken2 #


There is often interesting biological information ...blalba

We assume we have miniconda installed, and we have a miniconda environment called "kraken" in which kraken2 is installed.


First we activate the miniconda environment (named "kraken") in which we have installed the kraken software so we can use it:

```sh
conda activate kraken
```

Now we first have to download and build a database containing the *l*-mers and the taxonomic information. We are using the *kraken2-build --standard* command for that. This command will download and build the standard Kraken 2 database. Afterwards we can use that database to classify sequencing reads as belonging to specific taxa. You need to replace $DBNAME with the path to where you want to save the database and the name you want to give the database. You can also choose the number of threads you  want kraken2 to use (here we use 20, but you can change that, if you want). The database building will take a while.

```sh
kraken2-build --standard --threads 20 --db $DBNAME
```



