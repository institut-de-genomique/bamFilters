#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>
// lib hts
#include "htslib/sam.h"
// lib bamFilters
#include "getfa.h"

int nword = 16;
int nk = 3;
long int top_chrono;

void dust();
void dust_rect() {}
void set_dust_level(int value);
void getfafun(char *name, void (*fun)());
void usage();
void demarrer_chrono();
void stop_chrono();

int main(int argc, char** argv) {
	int verbose = 0;
	char *BAMfile = 0, *BAMoutfile = 0, *BAMoutfile_uniq = 0, *LISTREFNAME_OUTPUT = 0, *STATS_FILE = 0;
	FILE* fichierListRefSeq = NULL;
	FILE* handler_stats_file = NULL;
	int bamOutFile = 0;
// DEFAULT ARGUMENTS
	int identity_min = 0;
	int alignment_min_size = 0;
	int ratio = 100;
	int nomasked_len = 0;
	int intron_min_size = -1;
// COMPTEURS
	int nb_ali = 0;
	int nb_filter_selected = 0;
	int nb_filter_selected_uniq = 0;
// PARSE ARGUMENTS
	int c;
	// LIST OF TRETMENTS
	int identity_filter = 0;
	int percent_alignement_filter = 0;
	int percent_low_complexity_filter = 0;
	int number_high_complexity_filter = 0;
	int uniq_filter = 0;
	int get_ref_seq_name_only = 0;
	int produce_stats_file = 0;

	if(argc <= 1) {
		fprintf(stderr, "[bamFilters] number of args invalid\n");
		usage();
	}

	while ((c = getopt(argc, argv, "::b:o:i:a:r:n:s:u:y:z:hv")) != -1) {
		switch (c) {
			case 'b':
				BAMfile = optarg;
				break;
			case 'o':
				BAMoutfile = optarg;
				bamOutFile = 1;
				break;
			case 'i':
				identity_min = atoi(optarg);
				identity_filter = 1;
				break;
			case 'a':
				alignment_min_size = atoi(optarg);
				percent_alignement_filter = 1;
				break;
			case 'r':
				ratio = atoi(optarg);
				percent_low_complexity_filter = 1;
				break;
			case 'n':
				nomasked_len = atoi(optarg);
				number_high_complexity_filter = 1;
				break;
			case 's':
				intron_min_size = atoi(optarg);
				break;
			case 'u':
				BAMoutfile_uniq = optarg;
				uniq_filter = 1;
				break;
			case 'y':
				STATS_FILE = optarg;
				produce_stats_file = 1;
				break;
			case 'z':
				LISTREFNAME_OUTPUT = optarg;
				get_ref_seq_name_only = 1;
				break;
			case 'h':
				usage();
			case 'v':
				verbose = 1;
				break;
			default :
				fprintf(stderr, "[bamFilters] arguments problem !\n");
				usage();
		}
	}

// GESTION DES ARGUMENTS
	if(get_ref_seq_name_only && (uniq_filter || bamOutFile)){
		fprintf(stderr, "[bamFilters] arguments problem : -z no compatible with -o and -u \n");
		usage();
	}

// GESTION HANDLER INPUT/OUTPUT
	samFile *fp_in = NULL;
	samFile *fp_out = NULL;
	samFile *fp_out_uniq = NULL;
	bam1_t *b= NULL;
	sam_hdr_t *header = NULL;
	errno = 0;


	//ouverture du fichier BAM en INPUT
	if(0 == (fp_in = sam_open(BAMfile, "r"))){
		if (errno) perror("[bamFilters] BAMfile hts_open");
		else fprintf(stderr, "[bamFilters] BAMfile hts_open failed.\n");
		return errno;
	}

	// ne pas faire si on recupere que le nom des sequences de référence avec au moins un match
	if(!get_ref_seq_name_only){
		//ouverture du fichier BAM en OUTPUT
		if(bamOutFile){
			errno = 0;
			if(0 == (fp_out = sam_open(BAMoutfile, "wb"))){
				if (errno) perror("[bamFilters] BAMoutfile hts_open");
				else fprintf(stderr, "[bamFilters] BAMoutfile hts_open failed.\n");
				return errno;
			}
		}
		//ouverture du fichier BAM UNIQ en OUTPUT
		if(uniq_filter){
			errno = 0;
			if(0 == (fp_out_uniq = sam_open(BAMoutfile_uniq, "wb"))){
				if (errno) perror("[bamFilters] BAMoutfile_uniq hts_open");
				else fprintf(stderr, "[bamFilters] BAMoutfile_uniq hts_open failed.\n");
				return errno;
			}
		}
	}else{
		errno = 0;
	    fichierListRefSeq = fopen(LISTREFNAME_OUTPUT, "w");
	    if(fichierListRefSeq == NULL){
			if (errno) perror("[bamFilters] Ref file name failed to open");
			else fprintf(stderr, "[bamFilters] Ref file name failed to open.\n");
	    	return errno;
	    }
	}

	// ouverture du handler pour le fichier de stats
	if(produce_stats_file){
		errno = 0;
		handler_stats_file = fopen(STATS_FILE, "w");
		if(handler_stats_file == NULL){
			if (errno) perror("[bamFilters] Ref file name failed to open");
			else fprintf(stderr, "[bamFilters] Ref file name failed to open.\n");
			return errno;
		}
	}

// GESTION DU HEADER
	//lecture du header dans le BAM en INPUT
	if(verbose) fprintf(stderr, "Start : lecture du header dans le BAM\n");
	if(verbose) demarrer_chrono();

	errno = 0;
	if (0 == (header = sam_hdr_read(fp_in))){
		if (errno) perror("[bamFilters] sam_hdr_read");
		else fprintf(stderr, "[bamFilters] sam_hdr_read failed.\n");
		return errno;
	}
	if(verbose) fprintf(stderr, "End : lecture du header dans le BAM\n");
	if(verbose) stop_chrono();

	// ne pas faire si on recupere que le nom des sequences de référence avec au moins un match
	if(!get_ref_seq_name_only){
		//ecriture des headers dans les BAM en OUTPUT
		if(bamOutFile){
			if(verbose) fprintf(stderr, "Start : ecriture header dans bamOutput\n");
			if(verbose) demarrer_chrono();

			errno = 0;
			if(0 != (sam_hdr_write(fp_out, header))){
				if (errno) perror("[bamFilters] BAMoutfile sam_hdr_write");
				else fprintf(stderr, "[bamFilters] BAMoutfile sam_hdr_write failed.\n");
				return errno;
			}

			if(verbose) fprintf(stderr, "End : ecriture header dans bamOutput\n");
			if(verbose) stop_chrono();
		}
		if(uniq_filter){
			if(verbose) fprintf(stderr, "Start : ecriture header dans bamOutput_uniq\n");
			if(verbose) demarrer_chrono();

			errno = 0;
			if(0 != (sam_hdr_write(fp_out_uniq, header))){
				if (errno) perror("[bamFilters] BAMoutfile_uniq sam_hdr_write");
				else fprintf(stderr, "[bamFilters] BAMoutfile_uniq sam_hdr_write failed.\n");
				return errno;
			}

			if(verbose) fprintf(stderr, "End : ecriture header dans bamOutput_uniq\n");
			if(verbose) stop_chrono();
		}
	}
// DEBUT PARSE BAM FILE
	//initialisation du BAM
	if(verbose) fprintf(stderr, "Start : initialisation du BAM\n");
	if(verbose) demarrer_chrono();

	errno = 0;
	if(0 == (b = bam_init1())){
		if (errno) perror("[bamFilters] bam_init1");
		else fprintf(stderr, "[bamFilters] bam_init1 failed.\n");
		return errno;
	}
	if(verbose) fprintf(stderr, "End : initialisation du BAM\n");
	if(verbose) stop_chrono();
	if(verbose) fprintf(stderr, "Start : parcours du bam\n");
	if(verbose) demarrer_chrono();
	//printf("read_name\tread_size\tfrag_size\tsoft_clipping_size\thard_clipping_size\tpercent_ali\tnbr_match\tnbr_mis_match\tnbr_ins\tnbr_del\tsize_ali\tsum\tlow_comp\thigh_comp\tpercent_ide\tnbr_match_mismatch\tDE\n");
	//printf("read_name\tread_size\tpercent_ali\tpercent_ide\tnbr_low_compl_base\thigh_comp_nbr\tlow_comp_percent\n");
	int nbr_rm_id = 0;
	int nbr_rm_ali = 0;
	int nbr_rm_low = 0;
	int nbr_rm_high = 0;

	while(sam_read1(fp_in, header, b) >= 0) {
		if(b -> core.flag & BAM_FUNMAP){ // next record if unmapped read
			continue;
		}

		nb_ali++;
		if(nb_ali == 1000000){
			if(verbose) fprintf(stderr, "%d\n", nb_ali);
		}

		int match_mismatch=0;
		int match=0;
		int mismatch=0;
		int deletion=0;
		int insertion=0;
		int hard_clipping = 0;
		int soft_clipping=0;

		int taille_alignement;
		double percent_identity;
		int taille_fragment_lecture_aligne;
		int taille_lecture;
		double percent_alignement;

		// Count the mismatch number thanks SAM NM tag
		uint8_t *nm = bam_aux_get(b, "NM"); // default SAM format for edit distance of read (with star specific option --outSAMattributes and bwa default options)

		if (nm == NULL){
			fprintf(stderr, "[bamFilters] match without NM (edit distance) value. You need remove unmapped reads.\n");
			return 0;
            //nm = bam_aux_get(b, "nM"); // default STAR tag for edit distance of paire, not read.
		}

		int k;
		uint32_t *cigar = bam_get_cigar(b);

		// Count match + inser + del thanks the CIGAR code
		for(k=0;k< b->core.n_cigar;++k){
		  int cop =cigar[k] & BAM_CIGAR_MASK;//tag
		  int cl = cigar[k] >> BAM_CIGAR_SHIFT; //value
		  switch(cop){
			case BAM_CMATCH: match_mismatch+=cl;break; // nb match/mismatch
			case BAM_CINS: insertion+=cl;break; // nb insertions
			case BAM_CDEL: if (intron_min_size == -1 || cl < intron_min_size) deletion += cl; break; // nb deletions
			case BAM_CSOFT_CLIP: soft_clipping+=cl;break; // nb soft clipping
			case BAM_CHARD_CLIP: hard_clipping+=cl;break; // nb hard clipping
		  }
		}

		// taille alignement
		taille_alignement = match_mismatch + insertion + deletion;

		// nbr mismatch
		mismatch = nm[1];

		// nbr match
		match = match_mismatch - mismatch;

		percent_identity = ((double)match / taille_alignement)*100;

		taille_fragment_lecture_aligne = match + mismatch + insertion;

		taille_lecture = match + mismatch + insertion + soft_clipping + hard_clipping;

		percent_alignement = ((double)taille_fragment_lecture_aligne/taille_lecture)*100;

		if(produce_stats_file){
			char *readName = bam_get_qname(b);
			fprintf(handler_stats_file, "%50s\t",readName);
			fprintf(handler_stats_file, "%i\t",taille_lecture);
			fprintf(handler_stats_file, "%i\t",taille_fragment_lecture_aligne);
			fprintf(handler_stats_file, "%i\t",soft_clipping);
			fprintf(handler_stats_file, "%i\t",hard_clipping);
			fprintf(handler_stats_file, "%3.2f\t",percent_alignement);

			fprintf(handler_stats_file, "%i\t",match);
			fprintf(handler_stats_file, "%i\t",mismatch);
			fprintf(handler_stats_file, "%i\t",insertion);
			fprintf(handler_stats_file, "%i\t",deletion);
			fprintf(handler_stats_file, "%i\t",taille_alignement);
			fprintf(handler_stats_file, "%3.2f\t",percent_identity);
			fprintf(handler_stats_file, "\n");
		}


		// identity or/and % alignement filters
		if( identity_filter || percent_alignement_filter){
			// test % identity filter
			if(identity_filter && (percent_identity <= identity_min)){
				nbr_rm_id++;
				continue;
			}
			// test % alignement filter
			if(percent_alignement_filter && (percent_alignement <= alignment_min_size)){
				nbr_rm_ali++;
				continue;
			}
		}

		int bamseq_size;
		int sum;
		float frac;

		// Homofilter and/or Dustfilter filters
		if(percent_low_complexity_filter || number_high_complexity_filter){
			REGION *reg;
			char *bamseq = (char *) malloc(b->core.l_qseq+1);
			uint8_t *s = bam_get_seq(b);
			int n = 0;
			for(n=0;n<(b->core.l_qseq);n++) {
				int v = bam_seqi(s,n);
				bamseq[n] = seq_nt16_str[v]; // Table for converting a 4-bit encoded nucleotide to a letter.
			}
			bamseq[n] = 0;
			bamseq_size = strlen(bamseq);
			reg = dust_segs(bamseq_size, bamseq);
			int i = 0;
			sum = 0;
			for (i=0; reg[i].to != -1; i++) {
				sum += reg[i].to - reg[i].from;
			}
			frac = ((float)sum/bamseq_size)*100;

			free(bamseq);

			// test Homofilter : % low complexity filter
			if(percent_low_complexity_filter && (frac >= ratio)){
				nbr_rm_low++;
				continue;
			}

			// test Dustfilter : number of high-complexity bases
			if((number_high_complexity_filter && !((bamseq_size-sum) > nomasked_len))){
				nbr_rm_high++;
				continue;
			}
		}

		if(0){
			char *readName = bam_get_qname(b);
			fprintf(handler_stats_file, "%50s\t",readName);
			fprintf(handler_stats_file, "%i\t",taille_lecture);
			fprintf(handler_stats_file, "%i\t",taille_fragment_lecture_aligne);
			fprintf(handler_stats_file, "%i\t",soft_clipping);
			fprintf(handler_stats_file, "%i\t",hard_clipping);
			fprintf(handler_stats_file, "%3.2f\t",percent_alignement);

			fprintf(handler_stats_file, "%i\t",match);
			fprintf(handler_stats_file, "%i\t",mismatch);
			fprintf(handler_stats_file, "%i\t",insertion);
			fprintf(handler_stats_file, "%i\t",deletion);
			fprintf(handler_stats_file, "%i\t",taille_alignement);
			fprintf(handler_stats_file, "%3.2f\t",percent_identity);

			fprintf(handler_stats_file, "%i\t",sum); // nbr low complexity bases
			fprintf(handler_stats_file, "%i\t",(bamseq_size-sum)); // nbr high complexity base
			fprintf(handler_stats_file, "%3.2f\t",frac); // % low complexity

			fprintf(handler_stats_file, "%i\t",match_mismatch);
			fprintf(handler_stats_file, "%i\t",mismatch);
			fprintf(handler_stats_file, "%i\t",bamseq_size);
			fprintf(handler_stats_file, "\n");
		}

		// ne pas faire si on recupere que le nom des sequences de référence avec au moins un match
		if(!get_ref_seq_name_only){
			if(bamOutFile){
				sam_write1(fp_out, header, b);
				nb_filter_selected++;
			}

			// print uniq match in other bamFile
			if(uniq_filter){
				uint8_t *tag_uniq_bwa_aln = bam_aux_get(b, "X0");//uniquement avec bwa aln
				uint8_t *tag_XS = bam_aux_get(b, "XS");
				uint8_t *tag_AS = bam_aux_get(b, "AS");

				if (tag_uniq_bwa_aln != NULL && tag_uniq_bwa_aln[1] == 1) {
					sam_write1(fp_out_uniq, header, b);
					nb_filter_selected_uniq++;
				}else if(tag_XS != NULL && tag_AS != NULL){
					if(tag_AS[1] != tag_XS[1]){
						sam_write1(fp_out_uniq, header, b);
						nb_filter_selected_uniq++;
					}
				}

			}
		}else{
			char *seqid = header->target_name[b->core.tid]; // recuperation du nom de la sequence du match
			fprintf(fichierListRefSeq, "%s\n", seqid);

		}
	   // liberation de mémoire

	}
	if(verbose) fprintf(stderr, "[stats] rm for low identity : %i\n", nbr_rm_id);

	if(verbose) fprintf(stderr, "[stats] rm for low alignement : %i\n", nbr_rm_ali);

	if(verbose) fprintf(stderr, "[stats] rm for low complexity : %i\n", nbr_rm_low);

	if(verbose) fprintf(stderr, "[stats] rm for nbr base low complexity : %i\n", nbr_rm_high);

	if(verbose) fprintf(stderr, "[stats] Total reads : %i\n", nb_ali);

	if(verbose) fprintf(stderr, "[stats] Total rm reads : %i\n", (nbr_rm_id+nbr_rm_ali+nbr_rm_low+nbr_rm_high));

	if(verbose) fprintf(stderr, "End : parcours du bam\n");
	if(verbose) stop_chrono();
	// destruction des objets
	bam_destroy1(b);
	sam_hdr_destroy(header);
	// fermeture des handlers
	sam_close(fp_in);
	// ne pas faire si on recupere que le nom des sequences de référence avec au moins un match
	if(!get_ref_seq_name_only){
		if(bamOutFile){
			sam_close(fp_out);
		}

		if(uniq_filter){
			sam_close(fp_out_uniq);
		}

	}else{
		fclose(fichierListRefSeq);
	}

	if(produce_stats_file){
		fclose(handler_stats_file);
	}

	if(verbose) {
		fprintf(stderr, "[main_bamFilters] %d selected reads from %d input sequences.\n", nb_filter_selected, nb_ali);
		if (uniq_filter) {
			fprintf(stderr, "[main_bamFilters_uniq] %d selected unique reads from %d input sequences.\n", nb_filter_selected_uniq, nb_ali);
		}
	}
	return 0;
}

void demarrer_chrono() {
        top_chrono = clock();
}

void stop_chrono() {
        long int arret_chrono = clock();
        fprintf(stderr, "Le calcul a pris %f secondes.\n",
                (float)(arret_chrono - top_chrono) / CLOCKS_PER_SEC);
}

void usage() {
  fprintf(stderr, "--------------------------------------------------------------------------------------------\n");
  fprintf(stderr, "v-2019.01.22.0\n");
  fprintf(stderr, "compiled with dynamic htslib from samtools/1.10.2\n");
  fprintf(stderr, "Usage:   bamFilters -b <in.bam> -o <out.bma> [-u <out.uniq.bam>] [options] \n\n");
  fprintf(stderr, "Options: -b          FILE  input BAM file\n");
  fprintf(stderr, "         -o          FILE  output BAM file\n");
  fprintf(stderr, "         -i          INT   identity percent level, default is 0\n");
  fprintf(stderr, "         -a          INT   alignment percent level, default is 0\n");
  fprintf(stderr, "         -r          INT   filter out sequences which contain more than r%% of low-complexity bases, default is 100%%\n");
  fprintf(stderr, "         -n          INT   filter out sequences which contain less than n high-complexity bases, default is 0\n");
  fprintf(stderr, "         -s          INT   don't count deletions longer than s as parts of alignment (useful for spliced mapping)\n");
  fprintf(stderr, "         -u          FILE  output BAM file for uniq filter which select reads that mapped only at one position\n");
  fprintf(stderr, "         -z                bam files no genereted, only list of sequences names of reference with one or more match. \n                           "
		  "No compatible with -u and -o. print on stdout \n");
  fprintf(stderr, "         -y          FILE  produce stats files with percent id, ali for each alignement\n");
  fprintf(stderr, "         -v                verbose mode\n");
  fprintf(stderr, "         -h                help\n");
  fprintf(stderr, "         exemple pour tara : bamFilters -b ./in.bam -i 95 -a 80 -r 75 -n 30 -o ./out.bam -u ./out.uniq.bam \n");
  fprintf(stderr, "\n");
  fprintf(stderr, "--------------------------------------------------------------------------------------------\n");
  exit(1);
}
