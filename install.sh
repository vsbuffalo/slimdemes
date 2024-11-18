#!/bin/bash

# Installation script for SLiMDemes
# Downloads required Eidos files from vsbuffalo/slimdemes repository

# Set error handling
set -e

echo "Installing SLiMDemes..."

# Create directory if it doesn't exist
mkdir -p slimdemes

# Function to download a file
download_file() {
    local file=$1
    local url="https://raw.githubusercontent.com/vsbuffalo/slimdemes/master/eidos/${file}"

    echo "Downloading ${file}..."
    if ! curl -sS -o "slimdemes/${file}" "$url"; then
        echo "Error: Failed to download ${file}"
        exit 1
    fi
}

# Download required files
download_file "demes.eidos"
download_file "utilities.eidos"

# Check if files were downloaded successfully
if [ -f "slimdemes/demes.eidos" ] && [ -f "slimdemes/utilities.eidos" ]; then
    echo -e "\nslimdemes installed successfully!"
    echo -e "\nTo use slimdemes in your scripts, add these lines:"
    echo -e " source(\"slimdemes/demes.eidos\");"
    echo -e " source(\"slimdemes/utilities.eidos\");"
    echo -e "\nFor more information and documentation, visit:"
    echo "https://github.com/vsbuffalo/slimdemes"
else
    echo "Error: Installation failed"
    exit 1
fi
