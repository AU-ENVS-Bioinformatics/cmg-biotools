FROM bioperl/bioperl:release-1-7-2
RUN mkdir -p /usr/biotools

COPY biotools /usr/biotools

# Give execution permissions to everything in /usr/biotools and export it to the PATH
RUN chmod -R +x /usr/biotools/* && \
  echo "export PATH=$PATH:/usr/biotools" >> ~/.bashrc

# # Install Bio::SeqIO module perl
# #  Installing the dependencies failed: Module 'Graph::Directed' is not installed, Module 'XML::Parser::PerlSAX' is not installed, Module 'XML::DOM' is not installed, Module 'XML::Twig' is not installed, Module 'Test::Memory::Cycle' is not installed, Module 'XML::LibXML' is not installed, Module 'XML::LibXML::Reader' is not installed
## INstall GNU plot
RUN apt update && \
  apt upgrade -y && \
  apt install -y gnuplot lib32stdc++6
# Install PostScript::Simple.pm, Getopt::Long, 
RUN cpanm PostScript::Simple && \
  cpanm Getopt::Long && \
  cpanm Bio::SearchIO::blastxml

# Run bash in home directory
WORKDIR /root
CMD ["/bin/bash"]