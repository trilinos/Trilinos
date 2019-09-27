MAILCOMMAND="/usr/sbin/sendmail"
CONTENTTYPE="Content-Type: text/html"
RECIPIENTS=(
"wg-rol-trilinos@sandia.gov"
)
SUBJECT="Subject: ROL Test Summary $(date +%b-%d-%Y-at-%H:%M)\n\n"
TEXT="To avoid complacency, you must follow this link:\n\n"
LINK="https://testing.sandia.gov/cdash/index.php?project=Trilinos&date=now&filtercount=1&field1=site&compare1=63&value1=optimizer"

if [ "$1" -eq "1" ];then
  TEXT="There are test errors!!! Follow this link, pronto:\n\n"
fi

MESSAGE=${SUBJECT}${TEXT}${LINK}
echo -e ${MESSAGE} | ${MAILCOMMAND} ${RECIPIENTS[@]}

