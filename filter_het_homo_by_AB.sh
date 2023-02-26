#!/bin/bash

java -jar vcffilterjdk.jar \
-e 'return new VariantContextBuilder(variant).
genotypes(
    variant.
    getGenotypes().
    stream().
    map( G->{
     if(G.hasAD() && G.isHet()) {
            final int A[]= G.getAD();
            if( A.length>1 &&  A[0]>0 && A[1]/(double)A[0] > 0.249 && A[1]/(double)A[0] < 2.51 ) return G;
         }
        if(G.hasAD() && G.isHom()) {
            final int A[]= G.getAD();
         if( A[0] == 0 ) return G;
            if( A.length>1 &&  A[0]>0 && A[1]/(double)A[0] > 0.8999 ) return G;
            }
        if(G.hasAD() && G.isHomRef()) {
            final int A[]= G.getAD();
         if( A[1] == 0 ) return G;
            if( A.length>1 &&  A[1]>0 && A[1]/(double)A[0] < 0.11 ) return G;
            }
        return  GenotypeBuilder.createMissing(G.getSampleName(),2);
        }).
 collect(Collectors.toList())
    ).
make();' \
$1 \
| bgzip > $2
