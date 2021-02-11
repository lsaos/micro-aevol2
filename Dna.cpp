//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"

#include <cassert>

#ifdef OPTI_DNA_NO_MODULO

// Need to compute the maximum circular offset we will access
constexpr int MAX_DNA_MODULO(PROM_SIZE > TERM_SIZE ?
	(SHINE_DAL_SIZE + CODON_SIZE + SD_START_SPACER > PROM_SIZE ? SHINE_DAL_SIZE + CODON_SIZE + SD_START_SPACER : PROM_SIZE)
	: TERM_SIZE);

#endif // OPTI_DNA_NO_MODULO

Dna::Dna(int length, Threefry::Gen&& rng)
	: seq_(length
#ifdef OPTI_DNA_NO_MODULO
		+ MAX_DNA_MODULO
#endif // OPTI_DNA_NO_MODULO
	)
{
	// Generate a random genome
	for (int32_t i = 0; i < length; i++) {
#ifdef OPTI_DNA_BITSET
		seq_[i] = rng.random(NB_BASE) != 0;
#else // OPTI_DNA_BITSET
		seq_[i] = '0' + rng.random(NB_BASE);
#endif // OPTI_DNA_BITSET
	}

#ifdef OPTI_DNA_NO_MODULO
	std::copy(seq_.begin(), seq_.begin() + MAX_DNA_MODULO, seq_.begin() + length);
#endif // OPTI_DNA_NO_MODULO
}

void Dna::load(gzFile backup_file) {
	int dna_length;
	gzread(backup_file, &dna_length, sizeof(dna_length));

	std::vector<char> tmp_seq_buf(dna_length);
	char* const tmp_seq(tmp_seq_buf.data());

	gzread(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));

#ifdef OPTI_DNA_BITSET
	seq_ = SeqType(dna_length);
	for (int i = 0; i < dna_length; ++i) {
		if (tmp_seq[i] == '0') {
			seq_[i] = false;
		}
		else {
			seq_[i] = true;
		}
	}
#else // OPTI_DNA_BITSET
	seq_ = SeqType(tmp_seq, tmp_seq + dna_length);
#endif // OPTI_DNA_BITSET

#ifdef OPTI_DNA_NO_MODULO
	seq_.insert(seq_.end(), seq_.begin(), seq_.begin() + MAX_DNA_MODULO);
#endif // OPTI_DNA_NO_MODULO
}

void Dna::set(int pos, char c) {
#ifdef OPTI_DNA_BITSET
	if (c == '0') {
		seq_[pos] = false;
	}
	else {
		seq_[pos] = true;
	}
#else // OPTI_DNA_BITSET
	seq_[pos] = c;
#endif // OPTI_DNA_BITSET

#ifdef OPTI_DNA_NO_MODULO
	if (pos < MAX_DNA_MODULO) {
		seq_[(int)seq_.size() - MAX_DNA_MODULO + pos] = seq_[pos];
	}
#endif // OPTI_DNA_NO_MODULO
}

void Dna::do_switch(int pos) {
#ifdef OPTI_DNA_BITSET
	if (seq_[pos]) {
		seq_[pos] = false;
	}
	else {
		seq_[pos] = true;
	}
#else // OPTI_DNA_BITSET
	if (seq_[pos] == '0') {
		seq_[pos] = '1';
	}
	else {
		seq_[pos] = '0';
	}
#endif // OPTI_DNA_BITSET

#ifdef OPTI_DNA_NO_MODULO
	if (pos < MAX_DNA_MODULO) {
		seq_[(int)seq_.size() - MAX_DNA_MODULO + pos] = seq_[pos];
	}
#endif // OPTI_DNA_NO_MODULO
}

/**
 * Remove the DNA inbetween pos_1 and pos_2
 *
 * @param pos_1
 * @param pos_2
 */
void Dna::remove(int pos_1, int pos_2) {
	assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= length());

	seq_.erase(seq_.begin() + pos_1, seq_.begin() + pos_2);

#ifdef OPTI_DNA_NO_MODULO
	if (pos_1 < MAX_DNA_MODULO) {
		std::copy(seq_.begin(), seq_.begin() + MAX_DNA_MODULO, seq_.begin() + length());
	}
#endif // OPTI_DNA_NO_MODULO
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, SeqType seq) {
	// Insert sequence 'seq' at position 'pos'
	assert(pos >= 0 && pos < length());

	seq_.insert(seq_.begin() + pos, seq.begin(), seq.end());

#ifdef OPTI_DNA_NO_MODULO
	if (pos < MAX_DNA_MODULO) {
		std::copy(seq_.begin(), seq_.begin() + MAX_DNA_MODULO, seq_.begin() + length());
	}
#endif // OPTI_DNA_NO_MODULO
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, Dna* seq) {
	// Insert sequence 'seq' at position 'pos'
	assert(pos >= 0 && pos < length());

	seq_.insert(seq_.begin() + pos, seq->seq_.begin(), seq->seq_.end());

#ifdef OPTI_DNA_NO_MODULO
	if (pos < MAX_DNA_MODULO) {
		std::copy(seq_.begin(), seq_.begin() + MAX_DNA_MODULO, seq_.begin() + length());
	}
#endif // OPTI_DNA_NO_MODULO
}

void Dna::do_duplication(int pos_1, int pos_2, int pos_3) {
	// Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
	char* duplicate_segment = NULL;

	//int32_t seg_length;

	if (pos_1 < pos_2) {
		//
		//       pos_1         pos_2                   -> 0-
		//         |             |                   -       -
		// 0--------------------------------->      -         -
		//         ===============                  -         - pos_1
		//           tmp (copy)                      -       -
		//                                             -----      |
		//                                             pos_2    <-'
		//
		SeqType seq_dupl(seq_.begin() + pos_1, seq_.begin() + pos_2);

		insert(pos_3, seq_dupl);
	}
	else { // if (pos_1 >= pos_2)
	 // The segment to duplicate includes the origin of replication.
	 // The copying process will be done in two steps.
	 //
	 //                                            ,->
	 //    pos_2                 pos_1            |      -> 0-
	 //      |                     |                   -       - pos_2
	 // 0--------------------------------->     pos_1 -         -
	 // ======                     =======            -         -
	 //  tmp2                        tmp1              -       -
	 //                                                  -----
	 //
	 //
		SeqType seq_dupl(seq_.begin() + pos_1, seq_.end());
		seq_dupl.insert(seq_dupl.end(), seq_.begin(), seq_.begin() + pos_2);

		insert(pos_3, seq_dupl);
	}
}

int Dna::length() const {
	return (int)seq_.size()
#ifdef OPTI_DNA_NO_MODULO
		- MAX_DNA_MODULO
#endif // OPTI_DNA_NO_MODULO
		;
}

void Dna::save(gzFile backup_file) const {
	const int dna_length(length());

#ifdef OPTI_DNA_BITSET
	std::vector<char> buf(dna_length);
	for (int i = 0; i < dna_length; ++i) {
		if (seq_[i]) {
			buf[i] = '1';
		}
		else {
			buf[i] = '0';
		}
	}
#else // OPTI_DNA_BITSET
	const std::vector<char>& buf(seq_);
#endif // OPTI_DNA_BITSET

	gzwrite(backup_file, &dna_length, sizeof(dna_length));
	gzwrite(backup_file, buf.data(), dna_length * sizeof(buf[0]));
}

int Dna::promoter_at(int pos) const {
	int prom_dist[PROM_SIZE];

	for (int motif_id = 0; motif_id < PROM_SIZE; motif_id++) {
		int search_pos = pos + motif_id;

#ifndef OPTI_DNA_NO_MODULO
		if (search_pos >= (int)seq_.size()) {
			search_pos -= (int)seq_.size();
		}
#endif // OPTI_DNA_NO_MODULO

#ifdef OPTI_DNA_BITSET
		const bool checkVal(PROM_SEQ[motif_id] == '1');
#else // OPTI_DNA_BITSET
		const char checkVal(PROM_SEQ[motif_id]);
#endif // OPTI_DNA_BITSET

		// Searching for the promoter
		prom_dist[motif_id] = checkVal == seq_[search_pos] ? 0 : 1;
	}

	// Computing if a promoter exists at that position
	int dist_lead = prom_dist[0] +
		prom_dist[1] +
		prom_dist[2] +
		prom_dist[3] +
		prom_dist[4] +
		prom_dist[5] +
		prom_dist[6] +
		prom_dist[7] +
		prom_dist[8] +
		prom_dist[9] +
		prom_dist[10] +
		prom_dist[11] +
		prom_dist[12] +
		prom_dist[13] +
		prom_dist[14] +
		prom_dist[15] +
		prom_dist[16] +
		prom_dist[17] +
		prom_dist[18] +
		prom_dist[19] +
		prom_dist[20] +
		prom_dist[21];

	return dist_lead;
}

// Given a, b, c, d boolean variable and X random boolean variable,
// a terminator look like : a b c d X X !d !c !b !a
int Dna::terminator_at(int pos) const {
	int term_dist[TERM_STEM_SIZE];

	for (int motif_id = 0; motif_id < TERM_STEM_SIZE; motif_id++) {
		int right = pos + motif_id;
		int left = pos + (TERM_SIZE - 1) - motif_id;

#ifndef OPTI_DNA_NO_MODULO
		// loop back the dna inf needed
		if (right >= length()) {
			right -= length();
		}
		if (left >= length()) {
			left -= length();
		}
#endif // OPTI_DNA_NO_MODULO

		// Search for the terminators
		term_dist[motif_id] = seq_[right] != seq_[left] ? 1 : 0;
	}

	int dist_term_lead = term_dist[0] +
		term_dist[1] +
		term_dist[2] +
		term_dist[3];

	return dist_term_lead;
}

bool Dna::shine_dal_start(int pos) const {
	bool start = false;
	int t_pos, k_t;

	for (int k = 0; k < SHINE_DAL_SIZE + CODON_SIZE; k++) {
		k_t = k >= SHINE_DAL_SIZE ? k + SD_START_SPACER : k;
		t_pos = pos + k_t;

#ifndef OPTI_DNA_NO_MODULO
		if (t_pos >= (int)seq_.size()) {
			t_pos -= (int)seq_.size();
		}
#endif // OPTI_DNA_NO_MODULO

#ifdef OPTI_DNA_BITSET
		const char shineDalVal(SHINE_DAL_SEQ[k_t]);
		if (shineDalVal != '*' && seq_[t_pos] == (shineDalVal == '1')) {
#else // OPTI_DNA_BITSET
		if (seq_[t_pos] == SHINE_DAL_SEQ[k_t]) {
#endif // OPTI_DNA_BITSET
			start = true;
		}
		else {
			start = false;
			break;
		}
	}

	return start;
}

bool Dna::protein_stop(int pos) const {
	bool is_protein;
	int t_k;

	for (int k = 0; k < CODON_SIZE; k++) {
		t_k = pos + k;

#ifndef OPTI_DNA_NO_MODULO
		if (t_k >= (int)seq_.size()) {
			t_k -= (int)seq_.size();
		}
#endif // OPTI_DNA_NO_MODULO

#ifdef OPTI_DNA_BITSET
		const bool checkVal(PROTEIN_END[k] == '1');
#else // OPTI_DNA_BITSET
		const char checkVal(PROTEIN_END[k]);
#endif // OPTI_DNA_BITSET

		if (seq_[t_k] == checkVal) {
			is_protein = true;
		}
		else {
			is_protein = false;
			break;
		}
	}

	return is_protein;
}

int Dna::codon_at(int pos) const {
	int value = 0;

	int t_pos;

#ifdef OPTI_DNA_BITSET
	constexpr bool checkVal(true);
#else // OPTI_DNA_BITSET
	constexpr char checkVal('1');
#endif // OPTI_DNA_BITSET

	for (int i = 0; i < CODON_SIZE; i++) {
		t_pos = pos + i;

#ifndef OPTI_DNA_NO_MODULO
		if (t_pos >= (int)seq_.size()) {
			t_pos -= (int)seq_.size();
		}
#endif // OPTI_DNA_NO_MODULO

		if (seq_[t_pos] == checkVal) {
			value += 1 << (CODON_SIZE - i - 1);
		}
	}

	return value;
}
