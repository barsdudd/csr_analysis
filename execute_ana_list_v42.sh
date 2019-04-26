#!/bin/sh

#csr_list_v42 ana_list 8 64000 POT RC

for I in 57 59 62 67 70
do
    csr_list_v42 ana_list_$I 8 64000 POT RC &
done