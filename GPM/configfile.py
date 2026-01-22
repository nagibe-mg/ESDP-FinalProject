import gpm

username_pps = "johanna.kasischke@uni-bonn.de"  # likely your mail, all in lowercase
password_pps = "<your PPS password>"  # likely your mail, all in lowercase
username_earthdata = "johannakas"
password_earthdata = "esdp2026WiSe1"
base_dir = "/Users/johanna/Documents/Uni/MASTER/EDSP/finalproject/ESDP-FinalProject/GPM"  # On Windows: "C:\\Users\\<path\\to\\a\directory>\\GPM"
gpm.define_configs(
    #username_pps=username_pps,
    #password_pps=password_pps,
    username_earthdata=username_earthdata,
    password_earthdata=password_earthdata,
    base_dir=base_dir,
)