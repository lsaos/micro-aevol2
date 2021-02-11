//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <vector>
#include <zlib.h>

#include "Threefry.h"
#include "aevol_constants.h"

#ifdef OPTI_DNA_BITSET
static_assert(NB_BASE == 2, "Bitset optimization currently only supports NB_BASE=2");
#endif // OPTI_DNA_BITSET

class Dna {

public:
#ifdef OPTI_DNA_BITSET
	typedef std::vector<bool> SeqType; // std::vector<bool> is implemented as a bitset
#else // OPTI_DNA_BITSET
	typedef std::vector<char> SeqType;
#endif // OPTI_DNA_BITSET

	Dna() = default;

	Dna(const Dna& clone) = default;

	Dna(int length, Threefry::Gen&& rng);

	~Dna() = default;
	void load(gzFile backup_file);

	void set(int pos, char c);

	/// Remove the DNA inbetween pos_1 and pos_2
	void remove(int pos_1, int pos_2);

	/// Insert a sequence of a given length at a given position into the DNA of the Organism
	void insert(int pos, SeqType seq);

	/// Insert a sequence of a given length at a given position into the DNA of the Organism
	void insert(int pos, Dna* seq);

	void do_switch(int pos);

	void do_duplication(int pos_1, int pos_2, int pos_3);

public:
	int length() const;

	void save(gzFile backup_file) const;

	int promoter_at(int pos) const;

	int terminator_at(int pos) const;

	bool shine_dal_start(int pos) const;

	bool protein_stop(int pos) const;

	int codon_at(int pos) const;

private:
	SeqType seq_;
};
