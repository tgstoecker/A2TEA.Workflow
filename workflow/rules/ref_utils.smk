import os
import posixpath
import gzip
import shutil

from urllib.request import urlretrieve
from urllib.request import urlcleanup


extensions=[".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz", ".gtf.gz"]


def is_url(file):
    return (
        file.startswith("http")
        or file.startswith("ftp")
        or file.startswith("sftp")
    )


def is_valid_annotation_file(annotation):
    for ext in extensions:
        if annotation.endswith("{}".format(ext)):
            return True


def get_file_ext(annotation):
    for ext in extensions:
        if annotation.endswith("{}".format(ext)):
            return ext


def is_gzipped(file):
    return(
        file.endswith(".gz")
    )


def gunzip_annotation(annotation_file):
    with gzip.open(annotation_file, 'rb') as f_in:
        with open(annotation_file.replace(".gz", ""), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def handle_annotation(annotation, species):
    print(annotation)
    print(get_file_ext(annotation))
    urlcleanup()
    if is_valid_annotation_file(annotation) and not is_gzipped(annotation):
        if is_url(annotation):
            print("Downloading annotation")
            urlretrieve(annotation, 'resources/' + species + '.annotation')
        elif not is_url(annotation):
            print("Linking annotation")
            os.symlink(os.path.abspath(annotation), 'resources/' + species + '.annotation')
    if is_valid_annotation_file(annotation) and is_gzipped(annotation):
        if is_url(annotation):
            print("Downloading annotation")
            urlretrieve(annotation, 'resources/' + species + '.annotation.gz')
        elif not is_url(annotation):
            print("Linking annotation")
            os.symlink(os.path.abspath(annotation), 'resources/' + species + '.annotation.gz')
        print("Unpacking annotation")
        gunzip_annotation('resources/' + species + '.annotation.gz')
        print("Done!")


def handle_pep_fasta(fasta, species):
    if is_url(fasta):
        print("Downloading fasta")
        if is_gzipped(fasta):
            urlretrieve(fasta, 'resources/' + species + '.pep.fa.gz')
            print("Unpacking fasta")
            with gzip.open('resources/' + species + '.pep.fa.gz', 'rb') as f_in:
                with open('resources/' + species + '.pep.fa', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            urlretrieve(fasta, 'resources/' + species + '.pep.fa')
    elif not is_url(fasta):
        print("Gunzipping or linking fasta file already on file system to resources/ dir")
        if is_gzipped(fasta):
            with gzip.open(os.path.abspath(fasta), 'rb') as f_in:
                with open('resources/' + species + '.pep.fa', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            os.symlink(os.path.abspath(fasta), 'resources/' + species + '.pep.fa')


def handle_cdna_fasta(fasta, species):
    if is_url(fasta):
        print("Downloading fasta")
        if is_gzipped(fasta):
            urlretrieve(fasta, 'resources/' + species + '.cdna.fa.gz')
            print("Unpacking fasta")
            with gzip.open('resources/' + species + '.cdna.fa.gz', 'rb') as f_in:
                with open('resources/' + species + '.cdna.fa', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            urlretrieve(fasta, 'resources/' + species + '.cdna.fa')
    elif not is_url(fasta):
        print("Gunzipping or linking fasta file already on file system to resources/ dir")
        if is_gzipped(fasta):
            with gzip.open(os.path.abspath(fasta), 'rb') as f_in:
                with open('resources/' + species + '.cdna.fa', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            os.symlink(os.path.abspath(fasta), 'resources/' + species + '.cdna.fa')


def handle_gen_fasta(fasta, species):
    if is_url(fasta):
        print("Downloading fasta")
        if is_gzipped(fasta):
            urlretrieve(fasta, 'resources/' + species + '.gen.fa.gz')
            print("Unpacking fasta")
            with gzip.open('resources/' + species + '.gen.fa.gz', 'rb') as f_in:
                with open('resources/' + species + '.gen.fa', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            urlretrieve(fasta, 'resources/' + species + '.gen.fa')
    elif not is_url(fasta):
        print("Gunzipping or linking fasta file already on file system to resources/ dir")
        if is_gzipped(fasta):
            with gzip.open(os.path.abspath(fasta), 'rb') as f_in:
                with open('resources/' + species + '.gen.fa', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            os.symlink(os.path.abspath(fasta), 'resources/' + species + '.gen.fa')
