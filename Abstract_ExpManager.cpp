//
// Created by elturpin on 03/12/2020.
//

#ifdef _WIN32
#include <direct.h>
#define mkdir(f, r) _mkdir(f)
#define err(e, m) fprintf(stderr, "Error %s (code %d)", m, e)
#else // _WIN32
#include <err.h>
#endif // _WIN32

#include <sys/stat.h>
#include "Abstract_ExpManager.h"

/**
 * Create stats and backup directory
 */
void Abstract_ExpManager::create_directory() {
    // Backup
    int status = mkdir("backup", 0755);
    if (status == -1 && errno != EEXIST) {
        err(EXIT_FAILURE, "backup");
    }

    // Stats
    status = mkdir("stats", 0755);
    if (status == -1 && errno != EEXIST) {
        err(EXIT_FAILURE, "stats");
    }
}