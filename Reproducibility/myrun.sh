#!/bin/bash

## Test python installation
python -c 'import seaborn; import pandas' || {
  echo "Your Python installation is missing some required packages."
  echo "The following are required"
  echo ""
  echo "  seaborn"
  echo "  pandas"
  echo ""
  echo "If you are using the conda distribution, you can setup the exact same"
  echo "python environment used in the paper with the file `conda-environment.yml`"
  echo ""
  echo "Otherwise, just install the missing packages using your distribution's"
  echo "package manager or pip"
  exit 1
}


BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PATH=$PATH:$BASEDIR/../build/core:$BASEDIR/../scripts
NUM_RUNS=1

if [[ ! -f $BASEDIR/Data/dblp/dblp-vldb-publication.txt ]]
then
    echo "Unpacking DBLP dataset"
    bunzip2 $BASEDIR/Data/dblp/dblp-vldb-publication.txt.bz2
fi

##############################################################################
#
# Experiments for Figures 1, 2, 3, and 4
#
##############################################################################

function figures () {
  ## Run comparisons
  COMPARISON_RESULTS_DIR=comparison-results/movement2007lcc/km
  test -d $COMPARISON_RESULTS_DIR || mkdir $COMPARISON_RESULTS_DIR
  cd $COMPARISON_RESULTS_DIR
 
  #DATASET=$BASEDIR/Data/proteins/collins2007-lcc.txt
  #DATASET=$BASEDIR/Data/proteins/soc-gplus-lcc.txt
  #DATASET=$BASEDIR/Data/proteins/Email-Enron-lcc.txt
  #DATASET=$BASEDIR/Data/proteins/facebook.txt
  #DATASET=$BASEDIR/Data/dblp/dblp-vldb-publication.txt
  #DATASET=$BASEDIR/Data/proteins/gavin2006-lcc.txt
  #DATASET=$BASEDIR/Data/proteins/protein_network_data.txt
  #DATASET=$BASEDIR/Data/proteins/579138.protein.links.v10.5_graph_id.txt

  #500 1000 1500 2000 2500
  #20 40 60 80 100
  #1000 1500 2000 2500 3000
  #DATASET=$BASEDIR/Data/proteins/cilccbiomine.txt
  #DATASET=$BASEDIR/Data/proteins/Fruit-Fly.txt
  #DATASET=$BASEDIR/Data/proteins/graph.txt
  #DATASET=$BASEDIR/Data/proteins/gavin2006-lcc.txt
  #DATASET=$BASEDIR/Data/proteins/movement90_lcc.txt
  DATASET=$BASEDIR/Data/proteins/movement2007lcc.txt
  for TARGET in 600 700 800 900 1000
  do
      for RUN in $(seq $NUM_RUNS)
      do
        
	       #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 1 --random 1 --topk 0
         ugraph-Alo4 --graph $DATASET --target $TARGET --celf 1 --random 0 --topk 0
         #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 0 --random 0 -- topk 0
         #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 0 --random 0 --topk 1
    
    done
  done

  # DATASET=$BASEDIR/Data/proteins/graph_jaccard.txt
  # for TARGET in 20 40 60 80 100 
  # do
  #     for RUN in $(seq $NUM_RUNS)
  #     do
        
  #        #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 1 --random 1 --topk 0
  #        ugraph-Alo4 --graph $DATASET --target $TARGET --celf 1 --random 0 --topk 0
  #        #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 0 --random 0 -- topk 0
  #        #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 0 --random 0 --topk 1
    
  #   done
  # done

  # DATASET=$BASEDIR/Data/proteins/soc-gplus-lcc.txt
  # for TARGET in 1000 1500 2000 2500 3000
  # do
  #     for RUN in $(seq $NUM_RUNS)
  #     do
        
  #        #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 1 --random 1 --topk 0
  #        ugraph-Alo4 --graph $DATASET --target $TARGET --celf 1 --random 0 --topk 0
  #        #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 0 --random 0 -- topk 0
  #        #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 0 --random 0 --topk 1
    
  #   done
  # done

  # DATASET=$BASEDIR/Data/proteins/Email-Enron-lcc.txt
  # for TARGET in 1000 1500 2000 2500 3000
  # do
  #     for RUN in $(seq $NUM_RUNS)
  #     do
        
  #        #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 1 --random 1 --topk 0
  #        ugraph-Alo4 --graph $DATASET --target $TARGET --celf 1 --random 0 --topk 0
  #        #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 0 --random 0 -- topk 0
  #        #ugraph-Alo4 --graph $DATASET --target $TARGET --celf 0 --random 0 --topk 1
    
  #   done
  # done
}
  
## Run everything
figures

