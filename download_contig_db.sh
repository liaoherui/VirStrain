wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1oj-86Njz5mnY6djbhdv23a9r9OH5oqog' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1oj-86Njz5mnY6djbhdv23a9r9OH5oqog" -O VirStrain_contig_DB.tar.gz && rm -rf /tmp/cookies.txt  &&\

#https://drive.google.com/file/d/1oj-86Njz5mnY6djbhdv23a9r9OH5oqog/view?usp=drive_link

tar -zxvf VirStrain_contig_DB.tar.gz &&\
rm  VirStrain_contig_DB.tar.gz
