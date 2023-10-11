Bootstrap: docker
From: bioperl/bioperl:release-1-7-2

%post
    # Install necessary tools and configure environment
    apt update
    apt upgrade -y
    apt install -y gnuplot lib32stdc++6

    # Install Perl modules
    cpanm PostScript::Simple
    cpanm Getopt::Long
    cpanm Bio::SearchIO::blastxml

    # Copy Parallel directory to the appropriate location
    mkdir -p /usr/biotools
    cp -r biotools/* /usr/biotools
    chmod -R +x /usr/biotools/*
    # Install additional tools
    apt install -y ncbi-blast+ hmmer prodigal
%environment
    PATH=$PATH:/usr/biotools

%files
    biotools
