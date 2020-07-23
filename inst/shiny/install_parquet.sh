apt-get update
apt-get install curl wget gnupg
apt-get install -y -V apt-transport-https curl gnupg lsb-release
tee /etc/apt/sources.list.d/backports.list <<APT_LINE
deb http://deb.debian.org/debian $(lsb_release --codename --short)-backports main
APT_LINE
curl --output /usr/share/keyrings/apache-arrow-keyring.gpg https://dl.bintray.com/apache/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-keyring.gpg
tee /etc/apt/sources.list.d/apache-arrow.list <<APT_LINE
deb [arch=amd64 signed-by=/usr/share/keyrings/apache-arrow-keyring.gpg] https://dl.bintray.com/apache/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/ $(lsb_release --codename --short) main
deb-src [signed-by=/usr/share/keyrings/apache-arrow-keyring.gpg] https://dl.bintray.com/apache/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/ $(lsb_release --codename --short) main
APT_LINE
curl https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add -
tee /etc/apt/sources.list.d/llvm.list <<APT_LINE
deb http://apt.llvm.org/$(lsb_release --codename --short)/ llvm-toolchain-$(lsb_release --codename --short)-7 main
deb-src http://apt.llvm.org/$(lsb_release --codename --short)/ llvm-toolchain-$(lsb_release --codename --short)-7 main
APT_LINE
apt update
apt install -y -V libparquet-dev