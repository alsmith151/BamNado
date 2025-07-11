#!/bin/bash

# Pre-commit hooks setup script for BamNado
# This script installs pre-commit and sets up the hooks

set -e

echo "🔧 Setting up pre-commit hooks for BamNado..."

# Check if pre-commit is installed
if ! command -v pre-commit &> /dev/null; then
    echo "📦 Installing pre-commit..."

    # Try different installation methods
    if command -v pip &> /dev/null; then
        pip install pre-commit
    elif command -v pip3 &> /dev/null; then
        pip3 install pre-commit
    elif command -v brew &> /dev/null; then
        brew install pre-commit
    elif command -v conda &> /dev/null; then
        conda install -c conda-forge pre-commit
    else
        echo "❌ Could not install pre-commit. Please install it manually:"
        echo "   pip install pre-commit"
        echo "   or visit: https://pre-commit.com/#installation"
        exit 1
    fi
else
    echo "✅ pre-commit is already installed"
fi

# Install the git hook scripts
echo "🪝 Installing pre-commit hooks..."
pre-commit install

# Optionally install pre-push hooks for tests
echo "🚀 Installing pre-push hooks..."
pre-commit install --hook-type pre-push

# Run hooks on all files to check everything works
echo "🧪 Running hooks on all files to verify setup..."
pre-commit run --all-files || {
    echo "⚠️  Some hooks failed. This is normal on first run."
    echo "   The hooks have been installed and will run on future commits."
}

echo ""
echo "✅ Pre-commit hooks setup complete!"
echo ""
echo "📋 Available commands:"
echo "   pre-commit run --all-files    # Run all hooks on all files"
echo "   pre-commit run <hook-id>      # Run specific hook"
echo "   pre-commit autoupdate         # Update hook versions"
echo "   pre-commit uninstall          # Remove hooks"
echo ""
echo "🎯 The following hooks are now active:"
echo "   - cargo check (on Rust files)"
echo "   - cargo clippy (on Rust files)"
echo "   - cargo fmt --check (on Rust files)"
echo "   - cargo test --all-features (on push)"
echo "   - Basic file checks (trailing whitespace, YAML/TOML syntax, etc.)"
