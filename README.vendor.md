# Vendor Branch: vendor/molssi-cookiecutter

This branch serves as a **vendor branch** for importing and tracking the files from the [MolSSI cookiecutter](https://github.com/MolSSI/cookiecutter-cms) repository into [OpenAWSEM](https://github.com/cabb99/openawsem).
The main function of this branch is to allow changes to be tracked as parallel commits are done in the OpenAWSEM codebase and in the molssi-cookiecutter codebase.
It was created as an **orphan branch** to isolate the vendor files from the rest of the OpenAWSEM codebase:

## How to Update This Branch

To update the vendor files when the MolSSI cookiecutter makes changes:

1. **Check out this branch:**

   ```bash
   git checkout vendor/molssi-cookiecutter
   ```

2. **Run the cookiecutter**
    ``` bash
    $ cookiecutter gh:molssi/cookiecutter-cms
    ```
    ```
    You've downloaded /home/cb/.cookiecutters/cookiecutter-cms before. Is it okay to delete and re-download it? [y/n] (y): y
    [1/9] project_name (ProjectName): OpenAWSEM
    [2/9] repo_name (openawsem): 
    [3/9] first_module_name (openawsem): 
    [4/9] author_name (Your name (or your organization/company/team)): Carlos Bueno
    [5/9] author_email (Your email (or your organization/company/team)): carlos.bueno@rice.edu
    [6/9] description (A short description of the project (less than one line).): An implementation of the AWSEM (Associative memory, Water-mediated Structure, and Energy Model) coarse-grained protein forcefield designed for use with the OpenMM simulation toolkit.
    [7/9] Select open_source_license
        1 - MIT
        2 - BSD-3-Clause
        3 - LGPLv3
        4 - Not Open Source
        Choose from [1/2/3/4] (1): 1
    [8/9] Select dependency_source
        1 - Prefer conda-forge with pip fallback
        2 - Prefer default anaconda channel with pip fallback
        3 - Dependencies from pip only (no conda)
        Choose from [1/2/3] (1): 1
    [9/9] Select include_ReadTheDocs
        1 - y
        2 - n
        Choose from [1/2] (1): 1
    ```


2. **Remove the old vendor files:**

   ```bash
   git rm -rf .
   ```

3. **Copy the new versions of the cookiecutter files:**
   Replace the deleted files with their updated versions from the MolSSI cookiecutter output.

   ```bash
   rsync -av --progress COOKIECUTTER_FOLDER/openawsem/ . --exclude='.git'
   ```

4. **Revert this README**
    ```bash
    git checkout README.vendor.md
    ```

5. **Stage and commit the updates:**

   ```bash
   git add .
   git commit
   ```

---

## Merging into OpenAWSEM

Once the vendor branch is updated:

1. **Check out your integration branch** (e.g., `master`, or a feature branch in OpenAWSEM):

   ```bash
   git checkout master
   ```

2. **Merge or rebase the vendor branch:**

   ```bash
   git merge vendor/molssi-cookiecutter
   # or
   # git rebase vendor/molssi-cookiecutter
   ```

3. **Apply any further OpenAWSEM-specific adaptations** in follow-up commits.

---

## Notes

* This branch is **not meant for development**. Do not modify vendor files directly here.
* All adaptations or fixes should be done in the integration branch after the vendor merge.

---