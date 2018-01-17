#!/bin/sh
# runJobAgBiomeNextSeq.sh RUN=0RUN PLATE=0PLATE LANE=1|LANELIST=1,2,3,4 [STEP=mkRawNextSeq]
#    [ FilesDelivered=NO|YES CONT=[YES|NO] SAMPLESHEET=SampleSheet.csv STOP_AFTER=combine MINLEN=0 SUFFIX=_LANE7 REPLACE=NO|YES MVCP=cp ]
#    STEP="mkRawNextSeq|mkRaw|mkRaw384|mkRawMiSeq|replaceStep|mkTrim_Plate|mkTrim|mkAssemblies|copyEditConsensus|mkCLCAlignments|mkSam|clcDepthSampleSheet|runCLCReport"
# Use mkTrim_Plate!

# After the plate specific STEP=replaceStep, run: SampleSheet384to96.py; mkNewLane.sh L00#; runNewLane.sh
#
# After the plate specific STEP=runCLCReport there is a BREAK and typically STEP=combine is run
# STEP=combine is run ONLY once for all plates!
#
##    step="combine|combine_runJob_depth|combine_runJob_depth_plate|mkAssemblies_Velvet"
# Get list of plates that plate was run in with:
#    awk -F, '{if($2 ~/^A/){print $3}}' SampleSheet.csv  | sort -u
#
# WCC Sat Jan 11 12:43:21 CST 2014
#
# To update run 6:
# run replaceStep
# manually run on each Trim Plate: cd AgBiome_*_Trim/DATE_*_Plate##; idx2well.py ../../SampleSheet.csv
# rename each assembly directory: cd AgBiome_*_Assembly/DATE_*_Plate##; ~/bin/AgBmvAssemblyDir.sh Plate##
# run copyEditConsensus
# run mkCLCAlignments

echo ""

STARTTIME=`date`

PID=$$
echo "PID: $PID"

errFile=`basename $0 .sh`
errFile=${errFile}.err

CWD=`pwd`
echo "CWD: $CWD"

DATESTR=`basename ${CWD} | awk -F_ '{print $1}'`
#FLOWCELL=`basename ${CWD} | awk -F_ '{print $4}'`
FLOWCELL=`echo $CWD | awk -F'[/]' '{
   n=split($NF,a,"_");
   for(i=1;i<=n;i++){
      if (a[i] ~ /[AB]???????[xX][xX]/){
         print a[i];
         break;
      }
   }
}'`
echo "FLOWCELL: $FLOWCELL"

step=mkRawNextSeq
STOP_AFTER=""

FilesDelivered=NO

CONT="YES"
REPLACE=NO

OTHER_ARGS=""

SAMPLE_SHEET=SampleSheet.csv

MINLEN=0
SUFFIX=""

#MVCP=mv
MVCP=cp

LANE=""

echo ""
echo "cmd line: $0 $*"
echo ""

while [ $# -ge 1 ]
do
   name=`echo $1 | awk -F= '{print $1}'`
   value=`echo $1 | awk -F= '{print $2}'`

   if [ "$name" = "RUN" ]; then
      RUN=$value
   elif [ "$name" = "LANE" ]; then
      LANE=$value
   elif [ "$name" = "LANELIST" ]; then
      LANELIST=$value
   elif [ "$name" = "PLATE" ]; then
      PLATE=$value
   elif [ "$name" = "STEP" ]; then
      step=$value
   elif [ "$name" = "STOP_AFTER" ]; then
      STOP_AFTER=$value
   elif [ "$name" = "FilesDelivered" ]; then
      FilesDelivered=$value
   elif [ "$name" = "CONT" ]; then
      CONT=$value
   elif [ "$name" = "SAMPLESHEET" ]; then
      SAMPLE_SHEET=$value
   elif [ "$name" = "MINLEN" ]; then
      MINLEN=$value
   elif [ "$name" = "SUFFIX" ]; then
      SUFFIX=$value
   elif [ "$name" = "REPLACE" ]; then
      REPLACE=$value
   elif [ "$name" = "MVCP" ]; then
      # Default to moving Unaligned .fastq.gz to Raw. 
      # Support copying Unaligned .fastq.gz to Raw.
      MVCP=$value
   else
      OTHER_ARGS="$OTHER_ARGS $1"
      echo "added: $1 to OTHER_ARGS"
   fi
   
   shift
done

if [ -z "$LANE" ]; then
   if [ -n "$LANELIST" ]; then
      LANE=`echo $LANELIST | awk -F, '{print $1}'`
   else
      echo LANE and LANELIST unspecified
      exit 10
   fi
fi

# Flowcell_A seems to use Project_${FLOWCELL}, Flowcell_B: Project_Run${RUN}???
# See redifinition of FLOWCELL above
# that should work even if the run has been renamed
#DHOME=${CWD}/Unaligned_L00${LANE}/Project_Run${RUN}
#DHOME=${CWD}/Unaligned_L00${LANE}/Project_${FLOWCELL}
DHOME=${CWD}/Data/Intensities/BaseCalls
if [ "$step" != "demultiplex" ]; then
   if [ ! -d ${DHOME} ]; then
      if [ ! -d ${CWD}/Data/Intensities/BaseCalls ]; then
	  echo "Error: Could not find DHOME: $DHOME"
	  exit 1
      fi
   fi
fi

echo ""
echo "RUN: $RUN LANE: $LANE or LANELIST: $LANELIST PLATE: $PLATE step: $step DHOME: $DHOME"
JOBNAME=AgBiome_${DATESTR}_Run${RUN}
TITLE="AgBiome ${DATESTR} HiSeq Run ${RUN}: Plate ${PLATE}"

# NB: ${OUTDIR}/${JOBNAME}_Raw, _Trim 
OUTDIR=${CWD}

#SAV_index=${DHOME}/SAV_indexing.txt
#SAV_summary=${DHOME}/SAV_summary.txt
### AgBiome has been sending an Excel work book #paperwork=${CWD}/BCS119_scan.jpg
### AgBiome has not been registering on the website ??? #SRIDtxt=${CWD}/SRID1375.txt
###For HiSeq with 96 pooled samples, won't have individual pooled traces ##AGILENTLENtxt=${CWD}/BCS119_AgilentLen.txt

BT2_HOME=/usr/local/src/bowtie2-2.1.0

#### Does CLC cell have a trimmer? How good is it?
#FLOW_HOME=/usr/local/src/Partek/partek_flow

#CLCHOME=/opt/clc-assembly-cell-4.2.2-linux_64
#CLCHOME=/opt/clc-assembly-cell-4.3.0-linux_64
CLCHOME=/opt/clc-assembly-cell-5.0.2-linux_64
export PATH=$PATH:$CLCHOME

SAMHOME=/usr/local/src/samtools-0.1.19
export PATH=$PATH:$SAMHOME

VHOME=/usr/local/src/velvet_1.2.07
export PATH=$PATH:$VHOME

# Typically for Bayer use HASH_LEN=57
# but for AgBiome Velvet is a fall back assemblier so use HASH_LEN=21
#HASH_LEN=57
HASH_LEN=21
##EXP_COV=auto
EXP_COV=8
##COV_CUTOFF=2.448
COV_CUTOFF=auto

# The following has been replaced with bcl2fastq
#function demultiplex {
#   echo ""
#   echo "start demulitplex `date`"
#
#   for ((L=1; L<=8; L++))
#   do
#      echo "L: $L"
#      if [ ! -f SampleSheet_L00${L}.csv ]; then
#         if [ -f SampleSheet_L00${L}.csv ]; then
#            mysave SampleSheet_L00${L}.csv
#         fi
#         16to14SampleSheet.py ${SAMPLE_SHEET} l=${L} f=${FLOWCELL} > SampleSheet_L00${L}.csv
#      fi
#   
#      if [ ! -d ${CWD}/Unaligned_L00${L} ]; then
#         time /usr/local/bin/configureBclToFastq.pl \
#         --input-dir ${CWD}/Data/Intensities/BaseCalls \
#         --output-dir ${CWD}/Unaligned_L00${L} \
#         --sample-sheet ${CWD}/SampleSheet_L00${L}.csv \
#         --use-bases-mask 'Y*,I8,I8,Y*'  \
#         --flowcell-id $FLOWCELL \
#         --no-eamss \
#         --mismatches 1
#      fi
#   
#      cd Unaligned_L00${L}
#      nohup make -j 2 2>&1 > demultiplex.log &
#      cd ${CWD}
#   done
#
#   echo ""
#   echo "finish demultiplex `date`"
#}

function mkSam {
   echo ""
   echo "Start mkSam sample: ${s} `date`"

   adir=${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_Assembly
   if [ ! -d $adir ]; then
      echo "Could not find assembly dir: $adir"
      exit 1
   else
      ls -ltr $adir
   fi

   bn=`basename $ref _consensus${SUFFIX}.fasta`
   bn=${bn}${SUFFIX}

   # test -s file exists and has a size greater than zero
   if [ ! -s $adir/${bn}_sort_rmdup.bam ]; then

      if [ ! -s $adir/${bn}_sort.bam ]; then

         if [ ! -s $adir/${bn}.bam ]; then
	    echo "Convert cas to bam:"
            ${CLCHOME}/clc_cas_to_sam \
            -f 33 \
            -a $adir/${bn}.cas \
            -o $adir/${bn}.bam
         fi
         ls -l $adir/${bn}.bam

         echo "Sort the sam:"
         samtools sort $adir/${bn}.bam -T /tmp/${bn}_sorted -o $adir/${bn}_sort.bam
      fi
      if [ -s $adir/${bn}_sort.bam ]; then
         ls -l $adir/${bn}_sort.bam
	 if [ -s $adir/${bn}.bam ]; then
	     rm $adir/${bn}.bam
	 fi
      else
         echo "Error: Could not find .bam: $adir/${bn}_sort.bam"
         #exit 2
         return
      fi

      echo "Rmdup sam:"
      samtools rmdup $adir/${bn}_sort.bam $adir/${bn}_sort_rmdup.bam
   fi
   if [ -s $adir/${bn}_sort_rmdup.bam ]; then
      ls -l $adir/${bn}_sort_rmdup.bam
      if [ -s $adir/${bn}_sort.bam ]; then
	  rm $adir/${bn}_sort.bam
      fi
   else
      echo "Error: Could not find .bam: $adir/${bn}_sort_rmdup.bam"
      #exit 3
      return
   fi

   # This test won't be needed after one pass but doesn't hurt
   if [ -s $adir/${bn}_sort_rmdup_depth.txt ]; then
      rm $adir/${bn}_sort_rmdup_depth.txt
   fi

   if [ ! -s $adir/${bn}_sort_rmdup_stats.txt ]; then
      if [ -f $adir/${bn}_sort_rmdup.bam ]; then
         echo "start run depth2stats.py `date`"
         samtools depth $adir/${bn}_sort_rmdup.bam | \
         depth2stats.py name=${adir} | \
         tee -a summary.txt > $adir/${bn}_sort_rmdup_stats.txt
         echo "finish run depth2stats.py `date`"
      else
          echo "No: $adir/${bn}_sort_rmdup.bam found"
          echo "${bn}        0        0        0        0        0        0        0" > $adir/${bn}_sort_rmdup_stats.txt
      fi
   fi
   if [ -s $adir/${bn}_sort_rmdup_stats.txt ]; then
      ls -l $adir/${bn}_sort_rmdup_stats.txt
   else
      echo ""
      echo ""
      echo ""
      echo "Could not find .txt: $adir/${bn}_sort_rmdup_stats.txt"
      echo ""
      echo ""
      echo ""
      #exit 3
   fi

   echo ""
   echo "finish mkSam `date`"
}

function mkBTidx {
   echo ""
   echo "Start:  mkBTidx `date`"
   
   if [ ! -f ${reflist}.1.bt2 ]; then

      $BT2_HOME/bowtie2-build \
      ${reflist}.fasta \
      ${reflist}

      stat=$?
      if [ $stat != 0 ]; then
         echo Error in bowtie2-build $stat
         exit $stat
      fi
   fi

   echo ""
   echo "finish mkBTidx `date`"
}

function run_bowtie2 {
   echo ""
   echo "Start: run_bowtie2 `date`"

#   --very-senstive-local \
#   -k 2 \
#   -a \
#   --un-gz unpaired_reads_not_aligned.fastq.gz \
#   --al-gz unpaired_reads_aligned.fastq.gz \
#   --un-conc-gz paired_not_aligned_concordantly.fastq.gz \
#   --al-conc-gz paired_aligned_concordantly.fastq.gz \
#   --local \

   $BT2_HOME/bowtie2 \
   --time \
   --met-file metrics.txt \
   --threads 2 \
   --qc-filter \
   -x ${reflist} \
   -1 $trimReads1 \
   -2 $trimReads2 \
   -S ${s[$i]}_${i}_Assembly/${s[$i]}.sam
   stat=$?
   if [ $stat != 0 ];then
      echo "$stat Error bowtie2"
      echo "Trim can generate an unequal number of reads. For the purpose of getting TLEN, just continue"
#      exit $stat
   fi

   # -b for bam output
   # -h to include header
   # -o output BAM
   # -S input in SAM format
   # SAM file
   samtools view -b -h -o ${s[$i]}_${i}_Assembly/${s[$i]}.bam -S ${s[$i]}_${i}_Assembly/${s[$i]}.sam
   stat=$?
   if [ $stat != 0 ];then
      echo "$stat Error sam 2 bam"
      exit $stat
   fi

   if [ -s ${s[$i]}_${i}_Assembly/${s[$i]}.bam ]; then
      /bin/rm ${s[$i]}_${i}_Assembly/${s[$i]}.sam
   else
      echo "Could not find: ${s[$i]}_${i}_Assembly/${s[$i]}.bam"
      exit 10
   fi
   
   samtools sort ${s[$i]}_${i}_Assembly/${s[$i]}.bam -T /tmp/${s[$i]}_sorted -o ${s[$i]}_${i}_Assembly/${s[$i]}_sorted.bam
   stat=$?
   if [ $stat != 0 ];then
      echo "$stat Error sam sort"
      exit $stat
   fi

   if [ -s ${s[$i]}_${i}_Assembly/${s[$i]}_sorted.bam ]; then
      /bin/rm ${s[$i]}_${i}_Assembly/${s[$i]}.bam
   else
      echo "Could not find: ${s[$i]}_${i}_Assembly/${s[$i]}_sorted.bam"
      exit 11
   fi

   samtools rmdup ${s[$i]}_${i}_Assembly/${s[$i]}_sorted.bam ${s[$i]}_${i}_Assembly/${s[$i]}_sorted_rmdup.bam
   stat=$?
   if [ $stat != 0 ];then
      echo "$stat Error sam rmdup"
      exit $stat
   fi

   if [ -s ${s[$i]}_${i}_Assembly/${s[$i]}_sorted_rmdup.bam ]; then
      /bin/rm ${s[$i]}_${i}_Assembly/${s[$i]}_sorted.bam
   else
      echo "Could not find: ${s[$i]}_${i}_Assembly/${s[$i]}_sorted_rmdup.bam"
      exit 12
   fi

   # Find "good" pairs
   # To look at both reads for TLEN would be redundant so and only first read
   # http://samtools.sourceforge.net/samtools.shtml
   #samtools view -f pP1 -F uU2sfd ${s[$i]}_${i}_Assembly/${s[$i]}_sorted_rmdup.bam | SAMinsertHistogram.py

   requFlags=0x0043
   # 0x0001 = p = read is paired in sequencing
   # 0x0002 = P = read is mapped in a proper pair
   # 0x0040 = 1 = read is the first read in a pair

   skipFlags=0x078C
   # 0x0004 = u = query sequence itself is unmapped
   # 0x0008 = U = the mate is unmapped
   # 0x0080 = 2 = the read is the seecond read in a pair
   # 0x0100 = s = the alignment is not primary
   # 0x0200 = f = fails quality checks
   # 0x0400 = d = RCR or optical duplicate
   samtools view -f $reqFlages -F $skipFlags ${s[$i]}_${i}_Assembly/${s[$i]}_sorted_rmdup.bam | SAMinsertHistogram.py
   stat=$?
   if [ $stat != 0 ];then
      echo "$stat Error sam histogram"
      exit $stat
   fi

   echo ""
   echo "finish run_bowtie2 `date`"
}

function iReport {
   echo ""
   echo "start iReport `date`"
   IlluminaReport.py \
       title="${TITLE}" \
        runDir=${DHOME3} \
         index=${SAV_index} \
       summary=${SAV_summary} \
       paperwork=${paperwork} \
            SRID=${SRIDtxt} \
  AgilentLengths=${AGILENTLENtxt} \
        assemDir=${OUTDIR} \
          report=${OUTDIR}/report.htm

   echo ""
   echo "finish iReport `date`"
}

function mkRaw {
   echo ""
   echo "mkRaw: ${JOBNAME} Start: `date`"

   if [ ! -d ${JOBNAME}_Raw ]; then
      mkdir ${JOBNAME}_Raw
   fi
   if [ ! -d ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} ]; then
      mkdir ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}
   fi

   for s in `awk -F, -v flag=0 -v lane=${LANE} '{
      if (flag == 1) {
         if ($1 == lane) {
            print $3;
         }
      }
      if (flag == 2) {
         if ($2 == lane) {
            print $3;
         }
      }
      if ($0 ~ /^Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description/){
         flag=1;
      }
      if ($0 ~ /FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject/){
         flag=2;
      }
   }' ${SAMPLE_SHEET}`; 
   do 
      if [ ! -d ${DHOME}/Sample_${s} ]; then
         echo "mkRaw: $s not found"
      else
         for f in `ls ${DHOME}/Sample_${s}/${s}_*_L00${LANE}_R[12]_*.fastq.gz`
         do
            echo "f: $f"
            bname=`basename ${f} .fastq.gz`
            $MVCP ${f} ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname}_raw.fastq.gz
            stat=$?
            if [ $stat != 0 ]; then
               echo "Error moving ${f} ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname}_raw.fastq.gz"
               exit 1
            fi
            ls -l  ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname}_raw.fastq.gz
         done
      fi
   done

   echo ""
   echo "mkRaw: ${JOBNAME} End: `date`"
}

function mkRaw384 {
   echo ""
   echo "mkRaw384: ${JOBNAME} Start: `date`"

   if [ ! -d ${JOBNAME}_Raw ]; then
      mkdir ${JOBNAME}_Raw
   fi
   if [ ! -d ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} ]; then
      mkdir ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}
   fi

   FIRSTLANE=`echo $LANELIST | awk -F, '{print $1}'`
   for LANE in `echo $LANELIST | sed 's/,/ /g'`
   do

      for s_w in `awk -F, -v flag=0 -v lane=${LANE} -v plate=${PLATE} '{
         if (flag == 1) {
            if ($1 == lane && $4 ~ plate) {
               # Output the sample, well and index sequences
               # NB the space character in the printf at the end!!!!
               printf("%s_%s_%s_%s ",$3,$5,$7,$9);
            }
         }
         if (flag == 2) {
            # N.B. There is no PLATE to select from!!!!
            if ($2 == lane) {
               print $3;
            }
         }
         if ($0 ~ /^Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description/){
            # This is the MiSeq style sample sheet
            # I think this is the only case that works
            flag=1;
         }
         if ($0 ~ /FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject/){
            flag=2;
         }
      }' ${SAMPLE_SHEET}`; 
      do 
         s=`echo $s_w | awk -F_ '{print $1}'`
         w=`echo $s_w | awk -F_ '{print $2}'`
         i1=`echo $s_w | awk -F_ '{print $3}'`
         i2=`echo $s_w | awk -F_ '{print $4}'`
         echo "s_w: $s_w w: $w i1: $i1 i2: $i2"
         #ls ${CWD}/Unaligned_L00${LANE}/Project_${FLOWCELL}/Sample_${s}/${s}_${i1}-${i2}_L00${LANE}_R[12]_*.fastq.gz
         #exit 0

	 # Old demultiplex output dir:
	 fqDir=${CWD}/Unaligned_L00${LANE}/Project_${FLOWCELL}/Sample_${s}
	 # New demultiplex output dir:
	 fqDir=${CWD}/Unaligned_L00${LANE}
         
         if [ ! -d ${fqDir} ]; then
            echo "mkRaw384: $s for LANE: $LANE not found"
            continue
         else
	    # Old demultiplex output pat
            # for f in `ls ${fqDir}/${s}_${i1}-${i2}_L00${LANE}_R[12]_*.fastq.gz`
            # New demultiplex output pat
            for r in 1 2
            do
               for f in `ls ${fqDir}/${s}_S*_L00${LANE}_R${r}_*.fastq.gz`
               do
                  #bname=`basename ${f} .fastq.gz`
                  #bname2=`echo $bname | sed "s/L00[12345678]/L00${FIRSTLANE}/"`
                  bname2=${s}_${i1}-${i2}_L00${FIRSTLANE}_R${r}_001
                  echo "LANE: $LANE s: $s f: $f bname: $bname bname2: $bname2"
                  if [ -f  ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname2}_raw.fastq.gz ]; then
                     echo "found unexpected: ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname2}_raw.fastq.gz which we are appending to."
                     zcat ${f} ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname2}_raw.fastq.gz \
                     >> ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname2}_raw.fastq
                     rm ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname2}_raw.fastq.gz
                  else
                     zcat ${f} >> ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname2}_raw.fastq
                  fi
                  stat=$?
                  if [ $stat != 0 ]; then
                     echo "Error zcat ${f} ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname2}_raw.fastq"
                     exit 1
                  fi
                  ls -l  ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname2}_raw.fastq

               done
            done # for f in R1, R2
         fi
      done # for s in SampleSheet LANE
   done # for LANE in $LANELIST

   find ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} -name '*.fastq' -exec gzip {} \;

   # To allow the other functions that use: LANE to find the output:
   LANE=$FIRSTLANE

   echo ""
   echo "mkRaw384: ${JOBNAME} End: `date`"
}

function mkRawNextSeq {
   echo ""
   echo "mkRawNextSeq: ${JOBNAME} Start: `date`"

   if [ ! -d ${JOBNAME}_Raw ]; then
      mkdir ${JOBNAME}_Raw
   fi
   if [ ! -d ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} ]; then
      mkdir ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}
   fi

   for s_p_w in `awk -F, -v flag=0 -v lane=${LANE} -v plate=${PLATE} '{
      if (flag == 1) {
         if ($3 ~ plate) {
            # Output the sample, well and index sequences
            # NB the space character in the printf at the end!!!!
            #printf("%s_%s_%s_%s ",$2,$4,$6,$8);
            printf("%s_%s_%s_%s_%s ",$2,plate,$4,$6,$8);
         }
      }
      if (flag == 2) {
         # N.B. There is no PLATE to select from!!!!
         if ($2 == lane) {
            print $3;
         }
      }
      if (flag == 3) {
         if ($3 ~ plate) {
            # Output the sample, well and index sequences
            # NB the space character in the printf at the end!!!!
            #printf("%s_%s_%s_%s ",$2,$4,$6,$8);
            printf("%s_%s_%s_%s_%s ",$2,plate,$4,$6,$8);
         }
      }
      if ($0 ~ /^Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description/){
         flag=1;
      }
      if ($0 ~ /FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject/){
         flag=2;
      }
      if ($0 ~ /Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description/) {
         flag=3;
      }
      if ($0 ~ /Sample_ID,Sample_Name,Sample_Plate,Sample_96_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project/) {
         flag=3;
      }
   }' ${SAMPLE_SHEET}`; 
   do 
      s=`echo $s_p_w | awk -F_ '{print $1}'`
      p=`echo $s_p_w | awk -F_ '{print $2}'`
      w=`echo $s_p_w | awk -F_ '{print $3}'`
      i1=`echo $s_p_w | awk -F_ '{print $4}'`
      i2=`echo $s_p_w | awk -F_ '{print $5}'`
      echo "s_p_w: $s_p_w w: $w p: $p i1: $i1 i2: $i2"
      #ls ${CWD}/Unaligned_L00${LANE}/Project_${FLOWCELL}/Sample_${s}/${s}_${i1}-${i2}_L00${LANE}_R[12]_*.fastq.gz
      #exit 0


      fqDir=${CWD}/Data/Intensities/BaseCalls
      echo Demultiplex input dir: $fqDir
      
      if [ ! -d ${fqDir} ]; then
         echo "mkRawNextSeq: $s for Plate: $Plate not found"
         continue
      else
         for r in 1 2
         do
            for f in `ls ${fqDir}/${s}_S*_R${r}_*.fastq.gz`
            do
               #bname2=${s}_${i1}-${i2}_L00${FIRSTLANE}_R${r}_001
               bname2=${s}_${p}_${w}_L00${LANE}_R${r}_001
               echo "PLATE: $PLATE s: $s f: $f bname2: $bname2"
               $MVCP ${f} ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname2}_raw.fastq.gz
               stat=$?
               if [ $stat != 0 ]; then
                  echo "Error $MVCP ${f} ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname2}_raw.fastq.gz"
                  exit 1
               fi
               ls -l  ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname2}_raw.fastq.gz

            done

         done # for f in R1, R2
      fi

      #exit 0

   done # for s in SampleSheet LANE

echo ""
echo "mkRawNextSeq: ${JOBNAME} End: `date`"
}

function mkRawMiSeq {
   echo ""
   echo "mkRawMiSeq: ${JOBNAME} Start: `date`"

   if [ ! -d ${JOBNAME}_Raw ]; then
      mkdir ${JOBNAME}_Raw
   fi
   if [ ! -d ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} ]; then
      mkdir ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}
   fi

   for s in `awk -F, -v flag=0 -v sid=0 '{
      if (flag == 3) {
         sid+=1;
         printf("%s_%s ",$2,sid);
      }
      if ($0 ~ /Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description/) {
         flag=3;
      }
   }' ${SAMPLE_SHEET}` 
   do 
      sname=`echo $s | awk -F_ '{print $1}'`
      sid=`echo $s | awk -F_ '{print $2}'`
      echo "s: $s sname: $sname sid: $sid"
      if [ ! -f ${CWD}/Data/Intensities/BaseCalls/${sname}_S${sid}_L001_R1_001.fastq.gz ]; then
         echo "mkRawMiSeq: ${CWD}/Data/Intensities/BaseCalls/${sname}_S${sid}_L001_R1_001.fastq.gz not found"
         exit 1
      else
         for f in `ls ${CWD}/Data/Intensities/BaseCalls/${sname}_*_L00${LANE}_R[12]_*.fastq.gz`
         do
            echo "f: $f"
            bname=`basename ${f} .fastq.gz`
            cp ${f} ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname}_raw.fastq.gz
            stat=$?
            if [ $stat != 0 ]; then
               echo "Error moving ${f} ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname}_raw.fastq.gz"
               exit 1
            fi
            ls -l  ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname}_raw.fastq.gz
         done
      fi
   done

#      sname=`echo $s | awk -F_ '{print $1}'`
#      plate=`echo $s | awk -F_ '{print $2}'`
#      well=`echo $s | awk -F_ '{print $3}'`
#      sid=`echo $s | awk -F_ '{print $4}'`
#      echo "sname: $sname plate: $plate well: $well"
#      if [ -z "$plate" ]; then
#         plate=$PLATE
#         echo "Warning: SampleSheet.csv plate is blank. Using command line PLATE: $PLATE"
#      fi
#      if [ -z "$well" ]; then
#         echo "Well is blank check SampleSheet.csv expecting A1,B1,C1. Also check plate."
#         exit 15
#      fi
#      for r in R1 R2
#      do
#         if [ ! -f ${CWD}/Data/Intensities/BaseCalls/${sname}_S${sid}_${r}_001.fastq.gz ]; then
#            echo "mkRawMiSeq: ${CWD}/Data/Intensities/BaseCalls/${sname}_S${sid}_${r}_001.fastq.gz not found"
#            exit 10
#         else
#            if [ ! -f ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${sname}_${plate}_${well}_L001_${r}_001_raw.fastq.gz ]; then
#               cp ${CWD}/Data/Intensities/BaseCalls/${sname}_S${sid}_${r}_001.fastq.gz \
#                  ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${sname}_${plate}_${well}_L001_${r}_001_raw.fastq.gz
#               stat=$?
#               if [ $stat != 0 ]; then
#                  echo "Error copying ${CWD}/Data/Intensities/BaseCalls/${sname}_S${sid}_${r}_001.fastq.gz ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${sname}_${plate}_${well}_L001_${r}_001_raw.fastq.gz"
#                  exit 1
#               fi
#            fi
#            ls -l ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${sname}_${plate}_${well}_L001_${r}_001_raw.fastq.gz
#         fi
#      done
#   done


   echo ""
   echo "mkRawMiSeq: ${JOBNAME} End: `date`"
}

function run_trimbases {
   echo ""
   echo "Start: run_trimbases `date`"

   if [ -s $trimReads ]; then
      /bin/rm $trimReads
   fi

   #$FLOW_HOME/bin/trimbases-1.0/cutadapt \

   cutadapt \
   -q 20 \
   --minimum-length=${MINLEN} \
   $(</home/wcurtiss/Assembly/Illumina/Cutadapt/IlluminaAdaptorsNexTera.conf) \
   -o $trimReads \
   $rawReads 2>> testTrim_qual_adapt_5_3.err >> testTrim_qual_adapt_5_3.out
   stat=$?
   if [ $stat != 0 ];then
      echo "$stat Error cutadapt $rawReads1"
      exit $stat
   fi

   echo ""
   echo "finish run_trimbases `date`"
}

function mkTrim {
   echo ""
   echo "mkTrim: ${JOBNAME} Start: `date`"

   if [ ! -d ${JOBNAME}_Trim ]; then
      mkdir ${JOBNAME}_Trim
   fi

   if [ ! -d ${JOBNAME}_Trim/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} ]; then
      mkdir ${JOBNAME}_Trim/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}
   fi

   for s in `awk -F, -v flag=0 -v lane=${LANE} '{
      if (flag == 1) {
         if ($1 == lane) {
            print $3;
         }
      }
      if (flag == 2) {
         if ($2 == lane) {
            print $3;
         }
      }
      if (flag == 3) {
         if ($3 ~ plate) {
            # NB the space character in the printf at the end!!!!
            #printf("%s ",$2);
            printf("%s_%s_%s ",$2,$3,$4);
         }
      }
      if ($0 ~ /^Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description/){
         flag=1;
      }
      if ($0 ~ /FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject/){
         flag=2;
      }
      if ($0 ~ /Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description/) {
         flag=3;
      }
      if ($0 ~ /Sample_ID,Sample_Name,Sample_Plate,Sample_96_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project/) {
         flag=3;
      }
   }' ${SAMPLE_SHEET}`; 
   do 
         #for rawReads in `ls ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_*_L00${LANE}_R[12]_*.fastq.gz`
         for rawReads in `ls ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_L00${LANE}_R[12]_*.fastq.gz`
         do
            bname=`basename ${rawReads} _raw.fastq.gz`
	    trimReads=${JOBNAME}_Trim/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname}_trim.fastq.gz
	    if [ -s $trimReads ]; then
                if [ "$REPLACE" = "NO" ]; then
                   echo "Already found: $trimReads skipping"
		   continue
                else
                   echo "Already found: $trimReads replacing"
                   mysave $trimReads
                   rm $trimReads
                fi
	    fi
            run_trimbases

            ls -l $rawReads $trimReads
         done
#      fi
   done

   echo ""
   echo "mkTrim: ${JOBNAME} End: `date`"
}

function mkTrim_Plate {
   echo ""
   echo "mkTrim_Plate: ${JOBNAME} Start: `date`"

   if [ ! -d ${JOBNAME}_Trim ]; then
      mkdir ${JOBNAME}_Trim
   fi

   if [ ! -d ${JOBNAME}_Trim/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} ]; then
      mkdir ${JOBNAME}_Trim/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}
   fi

#####      if ($0 ~ /^Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description/){

   #for s_l in `awk -F, -v flag=0 -v plate=Ag-${PLATE} '{
   # Note use of function toupper() below to make it match Ag or AG
   for s in `awk -F, -v flag=0 -v plate=${PLATE} '{
      if (flag == 1) {
         if ($3 ~ plate) {
            # NB the space character in the printf at the end!!!!
            printf("%s ",$2);
         }
      }
      if ($0 ~ /Sample_ID/){
         flag=1;
      }
   }' ${SAMPLE_SHEET}`; 
   do 
         #s=`echo $s_l|awk -F_ '{print $1}'`

         for rawReads in `ls ${JOBNAME}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_*_L00${LANE}_R[12]_*.fastq.gz`
         do
            bname=`basename ${rawReads} _raw.fastq.gz`
	    trimReads=${JOBNAME}_Trim/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${bname}_trim.fastq.gz
	    if [ -s $trimReads ]; then
                if [ "$REPLACE" = "NO" ]; then
                   echo "Already found: $trimReads skipping"
		   continue
                else
                   echo "Alread found: $trimReads replacing"
                   mysave $trimReads
                   rm $trimReads
                fi
	    fi
            run_trimbases

            ls -l $rawReads $trimReads
         done
#      fi
   done

   echo ""
   echo "mkTrim_Plate: ${JOBNAME} End: `date`"
}

function mkAssemblies_Velvet {
   echo ""
   echo "start mkAssemblies Velvet `date`"

# Velvet assembly is a backup for CLC.
# Quit if the directory structure not already created because something is wrong.
   if [ ! -d ${JOBNAME}_Assembly ]; then
#      mkdir ${JOBNAME}_Assembly
      echo "No dir: ${JOBNAME}_Assembly found"
      exit 100
   fi

   if [ ! -d ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} ]; then
#      mkdir ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}
      echo "No dir:  ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} found"
      exit 101
   fi

#   for s in `selectSamples.py ${SAMPLE_SHEET} l=${LANE}`
   for s in `selectSamplesPlate.py ${SAMPLE_SHEET} p=${PLATE}`
   do 
      date

      interleaveStr=""
      for trimReads1 in `ls ${JOBNAME}_Trim/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_*L00${LANE}_R1_00?_trim.fastq.gz`
      do
          trimReads2=`echo $trimReads1 | sed 's/_R1_/_R2_/'`

          # Note that Velvet interleaveStr does not have the "-i" found in the CLC interleaveStr
          interleaveStr="$interleaveStr $trimReads1 $trimReads2"
   
          if [ -z "${trimReads1}" -o ! -f ${trimReads1} ]; then
             echo did not find valid $trimReads1 Continuing ....
             continue
          fi
   
          if [ -z "${trimReads2}" -o ! -f ${trimReads2} ]; then
             echo did not find valid $trimReads2 Continuing ....
             continue
          fi
      done

      adir=${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_Assembly
      if [ ! -d ${adir} ]; then
	  mkdir ${adir}
          #echo "No dir: ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_Assembly found"
          #exit 102
      fi

#      if [ -s ${adir}/${s}_${DATESTR}_consensus.fasta ]; then
#          echo already found ${adir}/${s}_${DATESTR}_consensus.fasta
#      else
          echo "Before Velvet Assembly ${s}: `date`"

          ${VHOME}/velveth \
          ${adir} \
          $HASH_LEN \
          -fastq.gz -shortPaired \
          -separate \
          $interleaveStr

#          $trimReads1 \
#          $trimReads2
    
          stat=$?

          echo "After velveth (step 1) assembly ${s}: `date` stat: $stat"
          if [ $stat != 0 ]; then
             echo error $stat in assembly: ${adir}/${s}_${DATESTR} | tee -a ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/Assembly.err
#             exit $stat
             continue
          fi
   
          ${VHOME}/velvetg \
          ${adir} \
          -exp_cov $EXP_COV \
          -cov_cutoff $COV_CUTOFF \
          -scaffolding yes
    
          stat=$?

          echo "After velvetg (step 2) assembly ${s}: `date` stat: $stat"
          if [ $stat != 0 ]; then
             echo error $stat in assembly: ${adir}/${s}_${DATESTR} | tee -a ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/Assembly.err
#             exit $stat
             continue
          fi

          inFa=${adir}/contigs.fa
          outFa=${adir}/AgBiome_${s}_${DATESTR}_consensus.fasta

          if [ -f $outFa ]; then
             mysave $outFa
             rm $outFa
          fi

          awk -F_ -v s=$s -v d=$DATESTR '{
             if ($0 ~ /^>/){
                printf(">%s_Contig%s_%s\n",s,$2,d)
             }else{
                print $0
             }
          }' $inFa | tee $outFa


#          if [ ! -s ${adir}/${s}_${DATESTR}_consensus.fasta ]; then
#             echo error missing or zero length ${adir}/${s}_${DATESTR}_consensus.fasta
#             ls -l ${adir}/${s}_${DATESTR}_consensus.fasta | tee -a ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/Assembly.err
#             continue
#          fi

#      fi

   done

   echo ""
   echo "finish mkAssemblies_Velvet `date`"
   
}

function mkAssemblies {
   echo ""
   echo "start mkAssemblies `date`"

   if [ ! -d ${JOBNAME}_Assembly ]; then
      mkdir ${JOBNAME}_Assembly
   fi

   if [ ! -d ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} ]; then
      mkdir ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}
   fi

   for s in `selectSamplesPlate.py ${SAMPLE_SHEET} p=${PLATE}`
   do 
      date

      interleaveStr=""
      for trimReads1 in `ls ${JOBNAME}_Trim/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_*L00${LANE}_R1_00?_trim.fastq.gz`
      do
          trimReads2=`echo $trimReads1 | sed 's/_R1_/_R2_/'`

          interleaveStr="$interleaveStr -i $trimReads1 $trimReads2"
   
          if [ -z "${trimReads1}" -o ! -f ${trimReads1} ]; then
             echo did not find valid $trimReads1 Continuing ....
             continue
          fi
   
          if [ -z "${trimReads2}" -o ! -f ${trimReads2} ]; then
             echo did not find valid $trimReads2 Continuing ....
             continue
          fi
      done

      adir=${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_Assembly
      if [ ! -d ${adir} ]; then
	  mkdir ${adir}
      fi

      if [ -s ${adir}/${s}_${DATESTR}_consensus.fasta ]; then
          echo already found ${adir}/${s}_${DATESTR}_consensus.fasta
      else
          echo "Before Assembly ${s}: `date`"

          # Argument order matters!!!! -p before -f -e
          #-q -i $trimReads1 $trimReads2 \

          ${CLCHOME}/clc_assembler \
          -v \
          --cpus 8 \
          -o ${adir}/${s}_${DATESTR}_consensus.fasta \
          -p fb ss 100 2000 \
          -q $interleaveStr \
          -f ${adir}/${s}_${DATESTR}.gff \
          -e ${adir}/${s}_${DATESTR}_dist.txt
    
          stat=$?

          echo "After assembly ${s}: `date` stat: $stat"
          if [ $stat != 0 ]; then
             echo error $stat in assembly: ${adir}/${s}_${DATESTR} | tee -a ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/Assembly.err
             #exit $stat
             continue
          fi

          if [ ! -s ${adir}/${s}_${DATESTR}_consensus.fasta ]; then
             echo error missing or zero length ${adir}/${s}_${DATESTR}_consensus.fasta
             ls -l ${adir}/${s}_${DATESTR}_consensus.fasta | tee -a ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/Assembly.err
             continue
          fi

      fi

   done

   echo ""
   echo "finish mkAssemblies `date`"
}

function mkAssembliesMetaGenomics {
   echo ""
   echo "start mkAssembliesMetaGenomics `date`"

   if [ ! -d ${JOBNAME}_Assembly ]; then
      mkdir ${JOBNAME}_Assembly
   fi

   if [ ! -d ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} ]; then
      mkdir ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}
   fi

   for s in `selectSamplesNamePlateWellMiSeq.py ${SAMPLE_SHEET} ${PLATE}`
   do 
      date

      interleaveStr=""
      for trimReads1 in `ls ${JOBNAME}_Trim/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_*L00${LANE}_R1_00?_trim.fastq.gz`
      do
          trimReads2=`echo $trimReads1 | sed 's/_R1_/_R2_/'`

          interleaveStr="$interleaveStr -i $trimReads1 $trimReads2"
   
          if [ -z "${trimReads1}" -o ! -f ${trimReads1} ]; then
             echo did not find valid $trimReads1 Continuing ....
             continue
          fi
   
          if [ -z "${trimReads2}" -o ! -f ${trimReads2} ]; then
             echo did not find valid $trimReads2 Continuing ....
             continue
          fi
      done

      adir=${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_Assembly
      if [ ! -d ${adir} ]; then
	  mkdir ${adir}
      fi

      if [ -s ${adir}/${s}_${DATESTR}_consensus.fasta ]; then
          echo already found ${adir}/${s}_${DATESTR}_consensus.fasta
      else
          echo "Before Assembly ${s}: `date`"

          # Argument order matters!!!! -p before -f -e
          #-q -i $trimReads1 $trimReads2 \

          ${CLCHOME}/clc_assembler \
          -v \
          --cpus 8 \
          -o ${adir}/${s}_${DATESTR}_consensus.fasta \
          -p fb ss 100 2000 \
          -q $interleaveStr \
          -f ${adir}/${s}_${DATESTR}.gff \
          -e ${adir}/${s}_${DATESTR}_dist.txt
    
          stat=$?

          echo "After assembly ${s}: `date` stat: $stat"
          if [ $stat != 0 ]; then
             echo error $stat in assembly: ${adir}/${s}_${DATESTR} | tee -a ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/Assembly.err
             #exit $stat
             continue
          fi

          if [ ! -s ${adir}/${s}_${DATESTR}_consensus.fasta ]; then
             echo error missing or zero length ${adir}/${s}_${DATESTR}_consensus.fasta
             ls -l ${adir}/${s}_${DATESTR}_consensus.fasta | tee -a ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/Assembly.err
             continue
          fi

      fi

   done

   echo ""
   echo "finish mkAssembliesMetaGenomics `date`"
}

function mkCLCAlignments {
   echo ""
   echo "start mkCLCAlignments `date`"

   if [ ! -d ${JOBNAME}_Assembly ]; then
      echo "missing Assembly dir: ${JOBNAME}_Assembly"
      exit 10
   fi

   if [ ! -d ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX} ]; then
      echo "missing Assembly Plate dir: ${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}"
      exit 11
   fi

   #for s in `selectSamplesPlateLanes.py ${SAMPLE_SHEET} p=${PLATE} l=${LANE}`
   for s in `selectSamplesPlate.py ${SAMPLE_SHEET} p=${PLATE}`
   do 
      date

      interleaveStr=""
      for trimReads1 in `ls ${JOBNAME}_Trim/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_*L00${LANE}_R1_00?_trim.fastq.gz`
      do
          trimReads2=`echo $trimReads1 | sed 's/_R1_/_R2_/'`

          interleaveStr="$interleaveStr -i $trimReads1 $trimReads2"

          if [ -z "${trimReads1}" -o ! -f ${trimReads1} ]; then
             echo did not find valid $trimReads1 Continuing ....
             continue
          fi
   
          if [ -z "${trimReads2}" -o ! -f ${trimReads2} ]; then
             echo did not find valid $trimReads2 Continuing ....
             continue
          fi
      done

      adir=${JOBNAME}_Assembly/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}/${s}_Assembly
      if [ ! -d ${adir} ]; then
	  echo "missing Assembly Plate Sample dir: ${adir}"
          #continue
	  exit 15
      fi
   
      ref=`ls ${adir}/AgBiome_${s}_${DATESTR}_consensus${SUFFIX}.fasta`
      if [ ! -f $ref ]; then
         echo "missing Assembly Plate Sample Consensus: $ref"
         continue
      fi

      bn=`basename $ref _consensus${SUFFIX}.fasta`
      bn=${bn}${SUFFIX}

#      for f in \
#${adir}/${s}_${DATESTR}.cas \
#${adir}/${s}_${DATESTR}_sort_rmdup.bam \
#${adir}/${s}_${DATESTR}_sort_rmdup_status.txt \
#${adir}/AgBiome_${s}_${DATESTR}.gff \
#${adir}/AgBiome_${s}_${DATESTR}.cas \
#${adir}/AgBiome_${s}_${DATESTR}_sort_rmdup.bam \
#${adir}/AgBiome_${s}_${DATESTR}_sort_rmdup_status.txt
#${adir}/AgBiome_${s}_${DATESTR}_sort_rmdup.bam \
#${adir}/AgBiome_${s}_${DATESTR}_sort_rmdup_status.txt
#      do
#          if [ -f ${f} ]; then
#              /bin/rm ${f}
#          fi
#      done

      if [ -s ${adir}/${bn}_sort_rmdup.bam ]; then
         echo "found result: ${adir}/${bn}_sort_rmdup.bam ... skipping"
         continue
      fi

      for f in \
${adir}/${bn}.cas \
${adir}/${bn}_sort_rmdup.bam \
${adir}/${bn}_sort_rmdup_status.txt \
${adir}/${bn}.gff
      do
          if [ -f ${f} ]; then
              /bin/rm ${f}
          fi
      done

      ${CLCHOME}/clc_mapper \
      --cpus 8 \
      -o ${adir}/${bn}.cas \
      -d ${ref} \
      -p fb ss 100 2000 \
      -q $interleaveStr

      #-d ${adir}/${s}_${DATESTR}_consensus.fasta \
      #-q -i $trimReads1 $trimReads2  

      stat=$?

      echo "After mapping ${s}: `date` stat: $stat"
      if [ $stat != 0 ]; then
         echo error $stat
         #exit $stat
         continue
      fi

      mkSam

   done
   echo ""
   echo "finish mkCLCAlignments `date`"
}

function mkBowtie2Alignments {
   echo ""
   echo "start mkBowtie2Alignments `date`"

   for ((i=${SSTRT}; i<=${SCNT}; i++))
   do
      date
   
      fq=`ls ${DHOME}/Data/Intensities/BaseCalls/*_S${i}_L001_R1_001.fastq.gz`; 
      sample=`basename $fq _S${i}_L001_R1_001.fastq.gz`; 
      echo ${i} ${sample} | tee -a sampleList.txt
      s[${i}]=${sample}
   
      rawReads1=${OUTDIR}/${JOBNAME}_Raw/${s[$i]}_S${i}_L001_R1_001_raw.fastq.gz
      rawReads2=${OUTDIR}/${JOBNAME}_Raw/${s[$i]}_S${i}_L001_R2_001_raw.fastq.gz
   
      trimReads1=${OUTDIR}/${JOBNAME}_Trim/${s[$i]}_S${i}_L001_R1_001_trim.fastq.gz
      trimReads2=${OUTDIR}/${JOBNAME}_Trim/${s[$i]}_S${i}_L001_R2_001_trim.fastq.gz
   
      if [ ! -d ${OUTDIR}/${s[$i]}_${i}_Assembly ]; then
          echo "Could not find: ${OUTDIR}/${s[$i]}_${i}_Assembly"
          exit 13
      fi
   
      reflist=${OUTDIR}/${s[$i]}_${i}_Assembly/contigs
      mkBTidx
   
      run_bowtie2
      if [ -s ${s[$i]}_${i}_Assembly/${s[$i]}_${i}_Assembly_consensus.1.bt2 ]; then
         rm ${s[$i]}_${i}_Assembly/${s[$i]}_${i}_Assembly_consensus.*.bt2
         rm ${s[$i]}_${i}_Assembly/${s[$i]}_sorted_rmdup.bam
      fi
      mkSam

   done
   echo ""
   echo "finish mkBowtie2Alginments `date`"
}
   
function combine_runJob_depth {
   cd ${CWD}

   if [ -f  AgBiome_${DATESTR}_Run${RUN}_Report/AgBiome_${DATESTR}_Run${RUN}.txt ]; then
      mysave  AgBiome_${DATESTR}_Run${RUN}_Report/AgBiome_${DATESTR}_Run${RUN}.txt
      /bin/rm  AgBiome_${DATESTR}_Run${RUN}_Report/AgBiome_${DATESTR}_Run${RUN}.txt
   fi

   hdrFlag=1
   for LANE in 1 2 3 4 5 6 7 8
   do
      if [ -f runJob_depth_L00${LANE}.txt ]; then
         if [ $hdrFlag == 1 ]; then
            head -1 runJob_depth_L00${LANE}.txt | \
            awk -F '[\t]' '{
               printf("Name\tRun Number\tSample Number\tLane Number\tPlate Number\tWell Number");
               for (i=2; i<=NF; i++) {
                  printf("\t%s",$i);
               }
               printf("\n");
            }' > \
            AgBiome_${DATESTR}_Run${RUN}_Report/AgBiome_${DATESTR}_Run${RUN}.txt
            hdrFlag=0
         fi
         #lc=`wc -l runJob_depth_L00${LANE}.txt | awk '{print $1}'`
         #let lc=${lc}-1
         #tail -${lc} runJob_depth_L00${LANE}.txt | \

         #tail -96 runJob_depth_L00${LANE}.txt | \
         cat runJob_depth_L00${LANE}.txt | \
         awk -F '[\t]' -v r=${RUN} -v l=${LANE} '{
            if ($0 !~ /^Name/) {
               n=split($1,a,"_");
               printf("%s\t%s\t%s\t%s\t%s\t%s",$1,r,NR,l,a[2],a[3]);
               for(j=2;j<=NF;j++){
                  printf("\t%s",$j)
               };
               printf("\n");
            }
         }' >> AgBiome_${DATESTR}_Run${RUN}_Report/AgBiome_${DATESTR}_Run${RUN}.txt
      fi
   done
}
   
function combine_runJob_depth_plate {
   # Purpose: to pull together runJob_depth_PLATE*.txt into a single report file.
   echo "start combine_runJob_depth_plate `date`"
   cd ${CWD}

   if [ -f  AgBiome_${DATESTR}_Run${RUN}_Report/AgBiome_${DATESTR}_Run${RUN}.txt ]; then
      mysave  AgBiome_${DATESTR}_Run${RUN}_Report/AgBiome_${DATESTR}_Run${RUN}.txt
      /bin/rm  AgBiome_${DATESTR}_Run${RUN}_Report/AgBiome_${DATESTR}_Run${RUN}.txt
   fi

   hdrFlag=1
   for PLATE in `awk -F, -v f=0 '{if(f==1){print $3};if($0 ~ /^Sample_ID/){f=1}}' $SAMPLE_SHEET  | sort -un`
   do
      #PLATE_LANE=`awk -F, -v p=$PLATE 'BEGIN{pp=sprintf("A[gG]-%d",p)}{if($4 ~ pp){print $1; exit}}' ${SAMPLE_SHEET}`
      #echo "PLATE: $PLATE PLATE_LANE: $PLATE_LANE"
      #echo "PLATE: $PLATE LANE: $LANE"
      echo "PLATE: $PLATE"
      if [ -f runJob_depth_PLATE${PLATE}.txt ]; then
         if [ $hdrFlag == 1 ]; then
            head -1 runJob_depth_PLATE${PLATE}.txt | \
            awk -F '[\t]' '{
               printf("Name\tRun Number\tSample Number\tLane Number\tPlate Number\tWell Number");
               for (i=2; i<=NF; i++) {
                  printf("\t%s",$i);
               }
               printf("\n");
            }' > \
            AgBiome_${DATESTR}_Run${RUN}_Report/AgBiome_${DATESTR}_Run${RUN}.txt
            hdrFlag=0
         fi

         cat runJob_depth_PLATE${PLATE}.txt | \
         awk -F '[\t]' -v r=${RUN} -v l=${LANE} '{
            if ($0 !~ /^Name/) {
               n=split($1,a,"_");
               printf("%s\t%s\t%s\t%s\t%s\t%s",$1,r,NR-1,l,a[2],a[3]);
               for(j=2;j<=NF;j++){
                  printf("\t%s",$j)
               };
               printf("\n");
            }
         }' >> AgBiome_${DATESTR}_Run${RUN}_Report/AgBiome_${DATESTR}_Run${RUN}.txt
      else
         echo "Warning: no runJob_depth_PLATE${PLATE}.txt found by combine_runJob_depth_plate"
      fi
   done
   echo "finish combine_runJob_depth_plate `date`"
}

function combine {
   echo ""
   echo "start combine `date`"
   #Run this once after all the lanes are finalized.
   #Create a directory that pulls all the consensus sequences into one folder:
   if [ ! -d AgBiome_${DATESTR}_Run${RUN}_Consensus ]; then
      mkdir AgBiome_${DATESTR}_Run${RUN}_Consensus
   fi

   find AgBiome_${DATESTR}_Run${RUN}_Assembly/ \
   -name 'AgBiome_*_consensus.fasta' \
   -exec cp {} AgBiome_${DATESTR}_Run${RUN}_Consensus \;

   if [ -f AgBiome_${DATESTR}_Run${RUN}_Consensus.tgz ]; then
      rm AgBiome_${DATESTR}_Run${RUN}_Consensus.tgz
   fi
   tar cvzf \
   AgBiome_${DATESTR}_Run${RUN}_Consensus.tgz  \
   AgBiome_${DATESTR}_Run${RUN}_Consensus

   if [ ! -d AgBiome_${DATESTR}_Run${RUN}_Consensus_Blast ]; then
      mkdir AgBiome_${DATESTR}_Run${RUN}_Consensus_Blast
   fi

   cd AgBiome_${DATESTR}_Run${RUN}_Consensus_Blast
   if [ -f AgBiome_${DATESTR}_Run${RUN}_Consensus_xml.log ]; then
      echo "Warning or Error: AgBiome_${DATASTR}_Run${RUN}_Consensus_xml.log is being replaced!!!!"
      mysave  AgBiome_${DATESTR}_Run${RUN}_Consensus_xml.log
      rm  AgBiome_${DATESTR}_Run${RUN}_Consensus_xml.log
   fi
   blastn_ncbi_16smicrobial_xml.sh \
      ${CWD}/AgBiome_${DATESTR}_Run${RUN}_Consensus

   if [ -f AgBiome_${DATESTR}_Run${RUN}_Consensus_table.txt ]; then
      mysave  AgBiome_${DATESTR}_Run${RUN}_Consensus_table.txt
      rm  AgBiome_${DATESTR}_Run${RUN}_Consensus_table.txt
   fi
   blast_xml.py \
      AgBiome_${DATESTR}_Run${RUN}_Consensus_xml.log > \
      AgBiome_${DATESTR}_Run${RUN}_Consensus_table.txt

   cd ${CWD}

   cd AgBiome_${DATESTR}_Run${RUN}_Report
   cat AgBiome_${DATESTR}_Run${RUN}_Plate*_Report.html > \
       AgBiome_${DATESTR}_Run${RUN}_Report.html

   
   # Pulled combine_runJob_depth out as a separate funciton 8/15/2014
   #combine_runJob_depth
   #combine_runJob_depth_plate

}

###################################################################

while test 1 -eq 1
do
#   if test "$step" = "demultiplex"; then
#      demultiplex
#      break
#      #step=mkRaw
#   elif test "$step" = "mkRaw"; then
   if test "$step" = "mkRaw"; then
      mkRaw
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "mkRaw" ]; then
         break
      else
          step="replaceStep"
      fi
   elif test "$step" = "mkRaw384"; then
      mkRaw384
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "mkRaw384" ]; then
         break
      else
          step="replaceStep"
      fi
   elif test "$step" = "mkRawNextSeq"; then
      mkRawNextSeq
      echo "run: /home/wcurtiss/bin/combineAgBiomeRaw.sh"
      break

#      if [ "$CONT" = "NO" ]; then
#         break
#      elif [ "$STOP_AFTER" = "mkRawNextSeq" ]; then
#         break
#      else
#          #step="replaceStep"
#          step="mkTrim_Plate"
#      fi
   elif test "$step" = "mkRawMiSeq"; then
      mkRawMiSeq
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "mkRawMiSeq" ]; then
         break
      else
          #step="replaceStep"
          step="mkTrim"
      fi
   elif test "$step" = "replaceStep"; then
      # Replace the Index sequences with the plate number and well number
      cd AgBiome_${DATESTR}_Run${RUN}_Raw/${DATESTR}_Run${RUN}_Plate${PLATE}${SUFFIX}
      #idx2well.py ../../SampleSheet.csv
      idx2well.py ../../${SAMPLE_SHEET}
      cd $CWD
      break

      #if [ "$CONT" = "NO" ]; then
      #   break
      #elif [ "$STOP_AFTER" = "replaceStep" ]; then
      #   break
      #else
      #   #step="mkTrim"
      #   step="mkTrim_Plate"
      #fi
   elif test "$step" = "mkTrim"; then
      mkTrim
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "mkTrim" ]; then
         break
      else
         step="mkAssemblies"
      fi
   elif test "$step" = "mkTrim_Plate"; then
      mkTrim_Plate
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "mkTrim_Plate" ]; then
         break
      else
         step="mkAssemblies"
      fi
   elif test "$step" = "mkAssemblies"; then
      mkAssemblies
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "mkAssemblies" ]; then
         break
      else
         step="copyEditConsensus"
      fi
   elif test "$step" = "mkAssembliesMetaGenomics"; then
      mkAssembliesMetaGenomics
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "mkAssembliesMetaGenomics" ]; then
         break
      else
         step="copyEditConsensus"
      fi
   elif test "$step" = "mkAssemblies_Velvet"; then
      mkAssemblies_Velvet
      break
#      if [ "$CONT" = "NO" ]; then
#         break
#      else
#         10/22/14 for RUN 20 LANE=8 AIM020248_154_E12
#         I did this manually because the contig number from Velvet would not be pulled out the same way.
#         step="copyEditConsensus"
#      fi
   elif test "$step" = "copyEditConsensus"; then
      echo "before `date` copyEditConsensus"
      copyEditConsensus.py DATESTR=$DATESTR RUN=$RUN PLATE=$PLATE LANE=$LANE SAMPLESHEET=$SAMPLE_SHEET SUFFIX=$SUFFIX
      stat=$?
      if [ $stat != 0 ]; then
         echo Error in copyEditConsensus.py $stat
         exit $stat
      fi
      echo "after `date` copyEditConsensus"
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "copyEditConsensus" ]; then
         break
      else
         step="mkCLCAlignments"
      fi
   elif test "$step" = "mkCLCAlignments"; then
      mkCLCAlignments
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "mkCLCAlignments" ]; then
         break
      else
         step="clcDepthSampleSheet"
      fi
   elif test "$step" = "mkSam"; then
      echo "Currently not implemented as stand alone. Needs to define s=sample in calling routine (clcDepthSampleSheet)"
      break
      mkSam
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "mkSam" ]; then
         break
      else
         step="clcDepthSampleSheet"
      fi
   elif test "$step" = "clcDepthSampleSheet"; then
      # This is where runJob_depth_L00#.txt are generated
      #clcDepthSampleSheetPlate.sh DATESTR=${DATESTR} RUN=${RUN} PLATE=${PLATE} LANE=${LANE} SUFFIX=${SUFFIX} SAMPLESHEET=${SAMPLE_SHEET}
      clcDepthSampleSheetPlate_v2.sh DATESTR=${DATESTR} RUN=${RUN} PLATE=${PLATE} LANE=${LANE} SUFFIX=${SUFFIX} SAMPLESHEET=${SAMPLE_SHEET}
      stat=$?
      if [ $stat != 0 ]; then
          echo "Error in: clcDepthSampleSheetPlate_v2.sh DATESTR=${DATESTR} RUN=${RUN} PLATE=${PLATE} LANE=${LANE} SUFFIX=${SUFFIX} SAMPLESHEET=${SAMPLE_SHEET}"
           exit $?
      fi
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "clcDepthSampleSheet" ]; then
         break
      else
         step="runCLCReport"
      fi
   elif test "$step" = "runCLCReport"; then
      #runCLCReport.sh $DATESTR $RUN $PLATE $LANE FilesDelivered=${FilesDelivered} ${SUFFIX}
      if [ -z "$LANELIST" ]; then
         #runCLCReport.sh $DATESTR $RUN $PLATE $LANE ${SUFFIX}
         runCLCReport_PLATE.sh $DATESTR $RUN $PLATE $LANE ${SUFFIX}
      else
         LANESAVE=$LANE
         for LANE in `echo $LANELIST | sed 's/,/ /g'`
         do
	    if [ -f runJob_depth_L00${LANE}.txt ]; then
               #runCLCReport.sh $DATESTR $RUN $PLATE $LANE ${SUFFIX}
               runCLCReport_PLATE.sh $DATESTR $RUN $PLATE $LANE ${SUFFIX}
            fi
         done
         LANE=$LANESAVE
      fi
      break
   elif test "$step" = "combine"; then
      combine
      if [ "$CONT" = "NO" ]; then
         break
      elif [ "$STOP_AFTER" = "combine" ]; then
         break
      else
         step="combine_runJob_depth_plate"
      fi
   elif test "$step" = "combine_runJob_depth"; then
      # I think this is where runJob_depth_L00#.txt are used
      combine_runJob_depth
      break
   elif test "$step" = "combine_runJob_depth_plate"; then
      # I think this is where runJob_depth_L00#.txt are used
      combine_runJob_depth_plate
      break
   elif test "$step" = "TEST"; then
      echo $*
      break
   else
      echo "unknown step: $step"
      break
   fi
done

FINISHTIME=`date`
echo "$FINISHTIME Finish"
echo "$STARTTIME Start"

