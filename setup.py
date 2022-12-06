import setuptools

with open("README.md","r") as fh:
    long_desc= fh.read()
    
setuptools.setup(
    name= "fit_to_class",
    version= "0.1",
    author= "mdebiaggi",
    author_email= "debiaggi@ianf.it",
    description= "Simple fitszille to classfits converter",
    long_description= long_desc,
    long_description_content_type= "text/markdown",
    url= "git@github.com:TheDebbio/fit_to_class.git",
    packages= ['commons','fit_to_class','fit_to_tests'],
    package_data={"fit_to_class": ["*.json"]},
    entry_points={
        "console_scripts":[
            "fit_to_class= fit_to_class.main:run"
            ]    
    }
)