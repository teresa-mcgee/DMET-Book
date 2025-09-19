from pathlib import Path
from openai import OpenAI
client = OpenAI()

# Paths
template_file = Path("/work/users/m/c/mcgeet/DMETpaper/DMET-results/example.qmd")
output_dir = Path("/work/users/m/c/mcgeet/DMETpaper/DMET-results")  # current directory, change if needed
biomarker_file = Path("/work/users/m/c/mcgeet/DMETpaper/DMET-results/biomarkernames_pyth.txt")

# Read template
content = template_file.read_text(encoding="utf-8")

# Read biomarker names from file (strip whitespace/quotes)
with biomarker_file.open("r", encoding="utf-8") as f:
    replacements = [line.strip().strip('"') for line in f if line.strip()]

# List of replacements
# replacements = ["Cyp2c39"]

# Create new files
for name in replacements:
    response = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": "You are a helpful scientific assistant."},
            {"role": "user", "content": f"Provide a very brief description of the mouse protein {name} collected from Uniprot, including gene name and aliases."}
        ]
    )
    gpt_reply = response.choices[0].message.content.strip()
    new_content = content.replace("{p}", name).replace("{GPT1}", gpt_reply)
    output_file = output_dir / f"{name}.qmd"
    output_file.write_text(new_content, encoding="utf-8")
    print(f"✅ Created {output_file}")
