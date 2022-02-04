#!/usr/bin/env bash
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1XYqr64tJec7VeDBD0Xc9cuUZqmawoty6' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1XYqr64tJec7VeDBD0Xc9cuUZqmawoty6" -O VirStrain_DB.tar.gz && rm -rf /tmp/cookies.txt  &&\

#https://drive.google.com/file/d/1XYqr64tJec7VeDBD0Xc9cuUZqmawoty6/view?usp=sharing

tar -zxvf VirStrain_DB.tar.gz &&\
rm  VirStrain_DB.tar.gz
