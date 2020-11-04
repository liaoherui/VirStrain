wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1KSVSOo5dc3i7KZZ2ReDaDNLt_Nf2DXw1' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1KSVSOo5dc3i7KZZ2ReDaDNLt_Nf2DXw1" -O VirStrain_DB.zip && rm -rf /tmp/cookies.txt  &&\

unzip VirStrain_DB.zip &&\
rm  VirStrain_DB.zip
