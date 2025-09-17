from pathlib import Path

# Paths
template_file = Path("/work/users/m/c/mcgeet/DMETpaper/DMET-results/example.qmd")
output_dir = Path("/work/users/m/c/mcgeet/DMETpaper/DMET-results")  # current directory, change if needed

# Read template
content = template_file.read_text(encoding="utf-8")

# List of replacements
replacements = ["Bsep", "Ent1"]

# Create new files
for name in replacements:
    new_content = content.replace("{p}", name)
    output_file = output_dir / f"{name}.qmd"
    output_file.write_text(new_content, encoding="utf-8")
    print(f"✅ Created {output_file}")
