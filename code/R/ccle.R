path.data <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/data/"


data <- read.table(paste(path.data, "CCLE_NP24.2009_Drug_data_2015.02.24.csv", 
                         sep=""), header=TRUE, stringsAsFactors=FALSE, sep=",")

data$CCLE.Cell.Line.Name
data2 <- data.frame(data, Tissue=unlist(lapply(strsplit(
  data$CCLE.Cell.Line.Name, "_"), function(l) {paste(l[-1], collapse=" ")})),
  stringsAsFactors=FALSE)


par(mar=c(19, 4, 4, 2) + 0.1)
barplot(table(data2$Tissue), las=2, ylab="Frequency", 
        main="Distribution of tissues")
par(mar=c(5, 4, 4, 2) + 0.1)

par(mar=c(7, 4, 4, 2) + 0.1)
barplot(table(data2$Compound), las=2, ylab="Frequency",
        main="Distribution of compounds")
par(mar=c(5, 4, 4, 2) + 0.1)

par(mar=c(5, 4, 4, 2) + 0.1)
barplot(sort(table(paste(data2$Compound, data2$Tissue)), decreasing=TRUE), 
        names.arg=FALSE, ylab="Freqency", xlab="Combinations",
        main="Distribution of compound and tissue combinations")
par(mar=c(5, 4, 4, 2) + 0.1)

hist(data2$IC50..uM., main="Distribution of IC50", xlab="IC50")


cell.lines <- unlist(strsplit("5637 	22Rv1 	786-O 	A-204 	A-253 	A2780 	A3/KAW 	A-375 	A4/Fuk 	A-498 	A549 	ABC-1 	AGS 	Alexander cells 	AN3-CA 	BGC-823 	BT-20 	BT-474 	C32 	CAL-12T 	CAL27 	Calu-1 	CHL-1 	CMK-86 	COLO-818 	COR-L105 	COR-L23 	DAN-G 	DB 	DBTRG-05MG 	DMS 79-CTG 	DOV13 	DU145 	EFE-184 	EFM-19 	EFO-27 	EN 	ESS-1 	EVSA-T 	FaDu 	G-401 	G-402 	GOS-3 	HCC1395 	HCC-366 	HCC-44 	HCC70 	HCC-78 	HCT-15 	HEC-1-A 	HEC-251 	HEC-265 	Hey-A8 	HLF 	Hs578T 	Hs766T 	Hs944.T 	HSC-2 	HT-1197 	HT-29 	HuH-1 	HuH-7 	HUP-T3 	HUP-T4 	HuTu 80 	IGR-39 	IGROV1 	IM95 	IMR-32 	J82 	JHH-1 	JHH-4 	JHOS-2 	JHUEM-1 	JHUEM-2 	KALS-1 	KATO III 	KG-1 	KM12 	KMBC-2 	KMRC-20 	KMS-11 	KMS-26 	KMS-27 	KMS-28BM 	KMS-34 	KNS-62 	KNS-81 	KP-1N 	KP-2 	KP-3 	KS-1 	KYSE-150 	KYSE-30 	KYSE-510 	L3.3 	LK-2 	LN-18 	LN-229 	LOU-NH91 	LoVo 	LS180 	M059K 	Malme3M 	MDA-MB-453 	MIA-PaCa-2 	MKN7 	MKN74 	MM1-S 	MOG-G-CCM 	MOLM-13 	MOLT-4 	MSTO-211H 	MV4-11 	NCI-H1299 	NCI-H1395 	NCI-H146-CTG 	NCI-H1563 	NCI-H1573 	NCI-H1650 	NCI-H1651 	NCI-H1793 	NCI-H1975 	NCI-H2052 	NCI-H2066 	NCI-H209-CTG 	NCI-H211 	NCI-H2122 	NCI-H2170 	NCI-H23 	NCI-H2347 	NCI-H2444 	NCI-H2452 	NCI-H441 	NCI-H446 	NCI-H522 	NCI-H524-CTG 	NCI-H661 	NCI-H810 	NCI-H838 	NCI-N87 	NIH:OVCAR-3 	NMC-G1 	NUGC-3 	NUGC-4 	ONS-76 	OV56 	OVCAR-8 	OVISE 	OVKATE 	OVSAHO 	OVTOKO 	Panc10.05 	PA-TU-8902 	PA-TU-8988T 	PC-3 	RERF-LC-Ad1 	RERF-LC-Ad2 	RERF-LC-MS 	RKO 	RMG-I 	RPMI-7951 	RT-112 	SBC-5 	SCaBER 	SCC-9 	SEM 	SJRH30 	SJSA-1 	SK-BR-3 	SK-CO-1 	SK-HEP-1 	SK-MEL-2 	SK-MEL-28 	SK-MEL-31 	SK-MEL-5 	SK-MES-1 	SK-N-BE 	SNB-19 	SNU-398 	SNU-423 	SW1088 	SW403 	SW48 	SW480 	SW579 	SW948 	T24 	T3M-4 	T-47D 	TCC-PAN2 	TE-11 	TE-9 	THP-1 	Toledo 	TOV-112D 	U-87MG 	U-937 	UACC-257 	VMRC-RCW 	WM-266-4 	YH13 	ZR-75-1 ", " \t"))
cell.lines2 <- data2[data2$Primary.Cell.Line.Name %in% cell.lines & data2$Compound=="Nutlin-3", c(2, 11)]
cell.lines3 <- cell.lines2[order(cell.lines2$Primary.Cell.Line.Name), ]
cell.lines3
