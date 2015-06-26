lune = window.lune = require 'lune'
console.log lune.phase()

game = new Phaser.Game window.innerWidth, window.innerHeight, Phaser.AUTO, '', {
  preload: ->
    game.load.image 'logo', 'luna.png'
  create: ->
    logo = game.add.sprite game.world.centerX, game.world.centerY, 'logo'
    logo.anchor.setTo 0.5, 0.5
}
