# Odin-TMB

Code to compute the TMB score from an input MAF file and assay (gene-panel) value.

## Version 0: DEPRECATED - Use with caution

This is the original version but it does not work properly in many cases. If the input MAF file has not been properly filtered to only have the events that contributed to the TMB count, then will get the wrong answer.

Ie the input MAF can only have events that are in the target area whose length was given as an argument.

