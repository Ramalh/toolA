# toolA
A tool proccesses simple biological data

Before using toolA, make sure "toolA" file is x(executible) permission. To check it write [ls -l] in the directoy where toolA is running. If it is not x mod, change the mod to x with [chmod +x toolA].

ToolA is the executible file, and it works in Linux command line and input file is fasta.

A command example: ./toolA --input fasta1.fasta --biomol_type DNA --GC_content_sub 5 --output output.txt (fasta file should be in the directory where command is running).

--input and -- biomol_type are only required options to fill.
