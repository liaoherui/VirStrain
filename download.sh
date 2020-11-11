wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=18IFeHuf85Qj5378QiFHGUJukxLV_RgYT' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=18IFeHuf85Qj5378QiFHGUJukxLV_RgYT" -O VirStrain_DB.tar.gz && rm -rf /tmp/cookies.txt  &&\

tar -zxvf VirStrain_DB.tar.gz &&\
rm  VirStrain_DB.tar.gz
