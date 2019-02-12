<?php

$hit_count = @file_get_contents('count-nix.txt');
$hit_count++;
@file_put_contents('count-nix.txt', $hit_count);

header('Location: http://iiirs.org/IIL/mtba/mtba-nix.zip'); // redirect to the real file to be downloaded

?>




