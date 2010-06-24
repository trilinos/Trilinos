#!/bin/tcsh

/bin/rm -f pure_output.*
set dir = `/bin/ls -d ch* hg*`
echo $dir
foreach d ($dir)
  echo $d
  (cd $d;  \
  grep -A 6 ABR: output/*.outerr >  ../pure_output.$d ; \
  grep -A 6 ABW: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 BRK: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 BSR: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 BSW: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 COR: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 FMM: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 FMR: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 FMW: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 FNH: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 FUM: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 IPR: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 IPW: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 MAF: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 MIU: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 MLK: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 MRE: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 MSE: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 NPR: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 NPW: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 PLK: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 SBR: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 SBW: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 SIG: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 SOF: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 UMC: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 UMR: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 WPF: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 WPM: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 WPN: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 WPR: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 WPW: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 WPX: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 ZPR: output/*.outerr >> ../pure_output.$d ; \
  grep -A 6 ZPW: output/*.outerr >> ../pure_output.$d ; \
  )
end
