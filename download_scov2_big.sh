wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1qAHjVADTiV3G00YekqystUXT2e7Ho2kq' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1qAHjVADTiV3G00YekqystUXT2e7Ho2kq" -O SCOV2_newBig.tar.gz && rm -rf /tmp/cookies.txt  &&\

#https://drive.google.com/file/d/1XYqr64tJec7VeDBD0Xc9cuUZqmawoty6/view?usp=sharing

tar -zxvf SCOV2_newBig.tar.gz &&\
rm  SCOV2_newBig.tar.gz
