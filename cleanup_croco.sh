#!/bin/bash
# cleanup_croco.sh — manage files in dated CROCO run folders
#
# Usage:
#   ./cleanup_croco.sh list    — show KEEP/DELETE status for every file
#   ./cleanup_croco.sh delete  — perform the actual cleanup
#
# Keep rules (per dated folder DD-MM-YYYY):
#   SCRATCH/croco_avg.nc
#   SCRATCH/croco_rst_*.nc
#   SCRATCH/d1_temp_salt_uv_z_all.nc
#   CROCO_FILES/**  (entire directory, EXCEPT croco_grd.nc)

set -euo pipefail

WORKSPACE="$(cd "$(dirname "$0")" && pwd)"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

usage() {
    echo "Usage: $(basename "$0") [list|delete]"
    echo "  list   — show KEEP/DELETE for every file in dated folders"
    echo "  delete — perform the actual cleanup"
    exit 1
}

# Return all dated folders (DD-MM-YYYY) sorted
get_dated_folders() {
    find "$WORKSPACE" -maxdepth 1 -type d \
        -regextype posix-extended \
        -regex '.*/[0-9]{2}-[0-9]{2}-[0-9]{4}' | sort
}

# is_kept <path-relative-to-dated-folder>
# Returns 0 (true) if the path should be kept, 1 (false) if it should be deleted
is_kept() {
    local rel="$1"

    # ---- SCRATCH files to keep ----
    [[ "$rel" == "SCRATCH/croco_avg.nc" ]]             && return 0
    [[ "$rel" == SCRATCH/croco_rst_*.nc ]]              && return 0
    [[ "$rel" == "SCRATCH/d1_temp_salt_uv_z_all.nc" ]]  && return 0
    [[ "$rel" == "SCRATCH/croco_forecast.in" ]]         && return 0

    # ---- CROCO_FILES: keep everything EXCEPT croco_grd.nc ----
    if [[ "$rel" == CROCO_FILES/* ]]; then
        [[ "$rel" == "CROCO_FILES/croco_grd.nc" ]] && return 1
        return 0
    fi

    # Everything else is deleted
    return 1
}

# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

do_list() {
    local keep_count=0 delete_count=0
    while IFS= read -r dated_dir; do
        local date_name
        date_name=$(basename "$dated_dir")
        echo ""
        echo "=== $date_name ==="
        while IFS= read -r file; do
            local rel="${file#"$dated_dir"/}"
            if is_kept "$rel"; then
                echo "  KEEP    $rel"
                (( keep_count++ )) || true
            else
                echo "  DELETE  $rel"
                (( delete_count++ )) || true
            fi
        done < <(find "$dated_dir" -type f 2>/dev/null | sort)
    done < <(get_dated_folders)
    echo ""
    echo "Total: $keep_count kept, $delete_count to delete"
}

do_delete() {
    local total=0
    while IFS= read -r dated_dir; do
        local date_name count=0
        date_name=$(basename "$dated_dir")
        echo ""
        echo "=== Cleaning $date_name ==="
        while IFS= read -r file; do
            local rel="${file#"$dated_dir"/}"
            if ! is_kept "$rel"; then
                echo "  Deleting $rel"
                rm -f "$file"
                (( count++ )) || true
            fi
        done < <(find "$dated_dir" -type f 2>/dev/null | sort)

        # Remove any directories that are now empty
        find "$dated_dir" -mindepth 1 -type d -empty -delete 2>/dev/null || true

        echo "  Deleted $count file(s)"
        (( total += count )) || true
    done < <(get_dated_folders)
    echo ""
    echo "Done. Total files deleted: $total"
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

case "${1:-}" in
    list)   do_list ;;
    delete) do_delete ;;
    *)      usage ;;
esac
