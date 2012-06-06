#!/bin/bash


infiles=$*


# Write infoline
echo "# D_Hb [um^2/s], D_Kni [um^2/s], Hb_mean [%EL], Hb_std_dev [%EL], Kni_mean [%EL], Kni_std_dev [%EL]"

# Initialize tempfile
echo -n > $0.temp

# Process input files
for infile in $infiles; do

# Grep internuclear spacing
l=$(cat $infile | grep "subv. size" | cut -d':' -f2)

# Grep diffusion constants
Dlist=$(cat $infile | grep "Diffusion" | cut -d':' -f2)
Drate_Hb=$( echo $Dlist | cut -d' ' -f29)
Drate_Kni=$(echo $Dlist | cut -d' ' -f30)

D_Hb=$( perl -e" print $l*$l*$Drate_Hb  / 4.0")
D_Kni=$(perl -e" print $l*$l*$Drate_Kni / 4.0")

cat $infile | awk -v D_Hb=$D_Hb -v D_Kni=$D_Kni '

	BEGIN{

		# Define the data file rows
		pos  = 4
		pHb  = 41	# Hb posterior boundary
		pKni = 74	# Kni anterior boundary

		# Init
		Hb_mean  = 0.0
		Hb_msqr  = 0.0

		Kni_mean = 0.0
		Kni_msqr = 0.0

		cc=""
	}

	{

		Hb_mean += $pos * $pHb
		Hb_msqr += $pos * $pos *$pHb

		Kni_mean += $pos * $pKni
		Kni_msqr += $pos * $pos *$pKni
	}

	END{

		Hb_var  = Hb_msqr  - Hb_mean*Hb_mean
		Kni_var = Kni_msqr - Kni_mean*Kni_mean

		sf = 100.0

		# Comment out the values for D=0
		if(D_Hb==0 || D_Kni==0)	cc="#"

		# Output is in [%EL]
		print cc, D_Hb, sf*Hb_mean, sf*sqrt(Hb_var), sf*Kni_mean, sf*sqrt(Kni_var)
	}

'	>> $0.temp

done

# Sort the values for ascending D_Hb
sort -n $0.temp

# Remove tempfile
rm $0.temp

