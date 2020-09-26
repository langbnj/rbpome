#!/users/gt/blang/bin/perl -w

# initialize

require('functions.inc.pl');
require('mysql.inc.pl');



# Download

# http://www.informatics.jax.org/downloads/reports/index.html#marker
# "MRK_List1.rpt (including withdrawn marker symbols)"
# "MRK_List2.rpt (excluding withdrawn marker symbols)" <<< This one (I don't want removed pseudogenes polluting the table)
# run("Download", "wget -O input/MRK_List1.txt 'http://www.informatics.jax.org/downloads/reports/MRK_List1.rpt'");
run("Download", "wget -O input/MRK_List2.txt 'http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt'");
# run("Unpack", "gunzip -f input/FILENAME");



# Doing this directly via NCBI homologene instead (the source)
# # Homolog mapping
# wget "http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt"
# wget "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
# wget "http://www.informatics.jax.org/downloads/reports/HGNC_homologene.rpt"

# vvvv Don't use this one, it's outdated. Doesn't contain this one: Eif1a MGI:95298 http://www.informatics.jax.org/marker/MGI:95298 According to this its human ortholog is EIF1AX but the table lists EIF1A. The other tables are fine.
# # wget "http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt"




# show directory
run("ls", "ls -lah input");

done();
