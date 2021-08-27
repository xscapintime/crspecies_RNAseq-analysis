#!/bin/bash

rsync -av --exclude={'*.txt','*.sh','*.pbs'} liyang@172.18.4.20:/home/liyang/data3/diapause/arrested_emb/streme/dngene/streme/ .
