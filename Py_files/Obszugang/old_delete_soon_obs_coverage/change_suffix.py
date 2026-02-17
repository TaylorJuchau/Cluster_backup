import os
from pathlib import Path

def change_file_suffixes(folder_path, old_suffix, new_suffix):
    """
    Changes the suffix of all files in a specified folder.

    Args:
        folder_path (str): The path to the directory containing the files.
        old_suffix (str): The current file extension (e.g., '.txt').
        new_suffix (str): The new file extension (e.g., '.log').
    """
    # Ensure suffixes start with a dot
    if not old_suffix.startswith('.'):
        old_suffix = '.' + old_suffix
    if not new_suffix.startswith('.'):
        new_suffix = '.' + new_suffix

    count = 0
    # Use pathlib to easily iterate and rename files
    for file_path in Path(folder_path).iterdir():
        # Check if it's a file and has the correct suffix
        if file_path.is_file() and file_path.suffix == old_suffix:
            # Construct the new filename
            new_file_path = file_path.with_suffix(new_suffix)
            try:
                # Rename the file
                file_path.rename(new_file_path)
                print(f"Renamed: {file_path.name} -> {new_file_path.name}")
                count += 1
            except OSError as e:
                print(f"Error renaming file {file_path.name}: {e}")

    print(f"\nOperation complete. Total files renamed: {count}")

# --- Example Usage ---

# 1. Define your parameters:
#    Replace 'path/to/your/folder' with the actual path to your folder.
#    On Windows, it might look like 'C:\\Users\\YourUser\\Desktop\\my_files'
#    On macOS/Linux, it might look like '/home/youruser/documents/my_files'
folder_to_process = 'data_output/'

#    Define the suffixes (include the dot for clarity)
current_suffix = '.npy'
target_suffix = '.pickle'

# 2. Run the function:
change_file_suffixes(folder_to_process, current_suffix, target_suffix)
