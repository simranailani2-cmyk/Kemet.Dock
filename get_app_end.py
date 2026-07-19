with open("app.py", "r") as f:
    lines = f.readlines()
print(f"Total lines: {len(lines)}")
for i, line in enumerate(lines[-30:]):
    print(f"{len(lines) - 30 + i + 1}: {line}", end="")
